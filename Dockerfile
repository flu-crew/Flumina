# Flumina container image
#
# Platform is pinned to linux/amd64 deliberately. gatk4, irma and lofreq have no
# linux-aarch64 builds on bioconda, so an arm64 image cannot be built at all.
# On Apple Silicon this runs under emulation (slower, correct); on any HPC node
# or cloud instance it runs natively.
#
# environment.yaml is the single source of truth for tool versions and is shared
# with the conda install path, so the container and a local conda env resolve to
# the same pinned versions.
#
#   docker build -t flumina:latest .
#   docker run --rm -v "$PWD":/data -w /data flumina:latest bwa 2>&1 | head -3

FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8-jammy

LABEL org.opencontainers.image.title="Flumina" \
      org.opencontainers.image.description="Variant calling pipeline for influenza A Illumina data" \
      org.opencontainers.image.source="https://github.com/flu-crew/Flumina" \
      org.opencontainers.image.licenses="MIT"

USER root
# procps supplies `ps`, which Nextflow shells out to for per-task resource metrics;
# without it every task logs a warning and reports no CPU/memory usage.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml \
 && micromamba clean --all --yes

# Nextflow launches processes with a plain `bash -ue`, which does not run the
# micromamba entrypoint, so the env must be on PATH unconditionally.
ENV PATH=/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# GATK spawns a JVM per task; without a cap it sizes the heap off the host's
# total RAM and gets OOM-killed inside a memory-limited container.
ENV JAVA_TOOL_OPTIONS="-XX:+UseSerialGC"

# Nextflow itself ships in the image so the container is self-sufficient: the
# user runs `apptainer run flumina_latest.sif ...` and everything, including the
# workflow engine, is already inside. Pinned rather than "latest" so rebuilding
# this Dockerfile cannot silently change the engine underneath the pipeline; the
# version must satisfy nextflow.config's manifest.nextflowVersion. curl and the
# JRE both come from the conda env above, so this needs nothing from apt.
# Installed after the conda layer so bumping it does not force that layer, which
# is the slow one, to rebuild.
USER root
# Pinned to match the Nextflow commonly available as a cluster module. This
# matters more than it looks: with `-p slurm` Nextflow runs on the HOST, so a
# cluster's version parses the config while the container's version runs the
# single-job path. When the two diverged (24.10.4 here, 26.04.3 on the cluster)
# every incompatibility surfaced only on the cluster, one at a time — the strict
# config parser introduced in 25.x rejects `def` functions, `if` statements and
# `for` loops, and stopped coercing `--flag false` to a boolean.
ARG NXF_VER=26.04.3
RUN NXF_VER=${NXF_VER} curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/nextflow \
 && chmod 755 /usr/local/bin/nextflow \
# The `nextflow` command is only a launcher: on first use it downloads the
# actual framework jar into NXF_HOME. Left to run time that download would
# repeat on every container start (NXF_HOME cannot persist in a read-only
# image) and would fail outright on an air-gapped cluster, so warm it here
# while the build still has network. 777 because the image is run by whatever
# uid the user happens to have — Apptainer never runs as root, and Docker is
# commonly run with -u — and Nextflow writes into NXF_HOME as it goes.
 && NXF_HOME=/opt/nextflow nextflow -version \
 && chmod -R 777 /opt/nextflow

# WFABC (Foll et al. 2015) is not packaged on conda or anywhere else, so it is
# built from the author's own source. Pinned to a commit rather than a branch so
# the binaries cannot change under a rebuild.
#
# No apt packages are needed: the conda env above already supplies make, a C++
# compiler, and the OpenMP runtime the Makefile's -fopenmp requires. Building
# with the conda toolchain also means the binaries link against the same
# libstdc++/libgomp that ship in this image.
USER root
ARG WFABC_COMMIT=c3d8896d91533fd0ca9f5a78fd6b5633dea8b90c
RUN cd /tmp \
 && curl -sSL "https://github.com/mfoll/WFABC/archive/${WFABC_COMMIT}.tar.gz" -o wfabc.tar.gz \
 && tar -xzf wfabc.tar.gz \
 && make -C "WFABC-${WFABC_COMMIT}/sources" CC=x86_64-conda-linux-gnu-g++ \
 && cp "WFABC-${WFABC_COMMIT}/sources/wfabc_1" \
       "WFABC-${WFABC_COMMIT}/sources/wfabc_2" /usr/local/bin/ \
 && chmod 755 /usr/local/bin/wfabc_1 /usr/local/bin/wfabc_2 \
 && rm -rf wfabc.tar.gz "WFABC-${WFABC_COMMIT}"

# IRMA from the CDC's own release. Bioconda stops at 1.0.3; the releases carry
# vendored blat, minimap2, pigz and parallel, so nothing else has to be added
# for it. Note IRMA's set_bin() prefers whatever is on PATH over its own
# vendored copies, which is why the parallel pin in environment.yaml still
# matters even though a good copy ships here.
ARG IRMA_VER=1.3.5
RUN cd /opt \
 && curl -sSL "https://github.com/CDCgov/irma/releases/download/v${IRMA_VER}/irma-v${IRMA_VER}-universal.zip" -o irma.zip \
 && unzip -q irma.zip \
 && rm irma.zip \
 && mv flu-amd irma \
 && chmod -R +x /opt/irma/IRMA /opt/irma/LABEL /opt/irma/IRMA_RES/third_party
ENV PATH=/opt/irma:$PATH

# SNPGenie, installed from the author's source rather than bioconda. The conda
# package depends on perl-list-util, which is built only against perl 5.22/5.26
# and therefore cannot coexist with irma's perl >=5.32 — the solver rejects the
# environment outright. The script itself needs nothing beyond perl core, so
# installing it directly sidesteps a conflict that has no other resolution.
# Pinned to a commit so a rebuild cannot silently change the analysis.
ARG SNPGENIE_COMMIT=d790569bc74622a64fbf5142e763581087bc7ea0
RUN cd /tmp \
 && curl -sSL "https://github.com/chasewnelson/SNPGenie/archive/${SNPGENIE_COMMIT}.tar.gz" -o snpgenie.tar.gz \
 && tar -xzf snpgenie.tar.gz \
# fasta2revcom.pl and gtf2revcom.pl are called by snpgenie.pl for minus-strand
# products, so they have to travel with it.
 && cp "SNPGenie-${SNPGENIE_COMMIT}/snpgenie.pl" \
       "SNPGenie-${SNPGENIE_COMMIT}/fasta2revcom.pl" \
       "SNPGenie-${SNPGENIE_COMMIT}/gtf2revcom.pl" /usr/local/bin/ \
 && chmod 755 /usr/local/bin/snpgenie.pl \
              /usr/local/bin/fasta2revcom.pl \
              /usr/local/bin/gtf2revcom.pl \
 && rm -rf snpgenie.tar.gz "SNPGenie-${SNPGENIE_COMMIT}"

# The pipeline itself. Copied last because it changes far more often than the
# tool environment above, so edits to it do not invalidate the conda layer.
COPY main.nf nextflow.config /opt/flumina/
# conf/ holds the settings shared by the scheduler profiles; nextflow.config
# includeConfig's it, so the image is unusable without it.
COPY conf /opt/flumina/conf
COPY reference.fa curated_database.csv /opt/flumina/
COPY example_file_rename.csv example_metadata.csv config.cfg irma_config.sh /opt/flumina/
COPY job_script_example_config.sh job_script_example_arguments.sh /opt/flumina/
COPY Scripts /opt/flumina/Scripts
RUN chmod -R +x /opt/flumina/Scripts

# Stamp when this image was built. The version string alone cannot distinguish
# two builds of the same version, so a stale cached .sif reports exactly what a
# fresh one does — which is a genuinely expensive thing to debug remotely. This
# sits after the COPY layers above, so any change to the pipeline invalidates it
# and the stamp is rewritten.
RUN date -u '+%Y-%m-%d %H:%M UTC' > /opt/flumina/BUILD_DATE
USER $MAMBA_USER

# FLUMINA_CONTAINER tells the launcher it is already inside the image, so it
# defaults to the `standard` profile and runs tasks in this container rather
# than trying to start a nested one.
#
# NXF_HOME points at the pre-warmed framework installed above. NXF_OFFLINE stops
# Nextflow reaching out for plugin and version checks at run time, which would
# otherwise fail on an air-gapped cluster and make every run depend on the
# network being up.
ENV FLUMINA_HOME=/opt/flumina \
    FLUMINA_CONTAINER=1 \
    NXF_HOME=/opt/nextflow \
    NXF_OFFLINE=true \
    PATH=/opt/flumina/Scripts:$PATH

WORKDIR /data

ENTRYPOINT ["flumina"]
