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
ARG NXF_VER=24.10.4
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

# The pipeline itself. Copied last because it changes far more often than the
# tool environment above, so edits to it do not invalidate the conda layer.
COPY main.nf nextflow.config /opt/flumina/
COPY reference.fa curated_database.csv /opt/flumina/
COPY example_file_rename.csv example_metadata.csv config.cfg irma_config.sh /opt/flumina/
COPY Scripts /opt/flumina/Scripts
RUN chmod -R +x /opt/flumina/Scripts
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
