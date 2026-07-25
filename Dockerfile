# Flumina container image
#
# Platform is pinned to linux/amd64 deliberately. gatk4, irma and lofreq have no
# linux-aarch64 builds on bioconda, so an arm64 image cannot be built at all.
# On Apple Silicon this runs under emulation (slower, correct); on any HPC node
# or cloud instance it runs natively.
#
# environment.yaml is the single source of truth for tool versions and is shared
# with the conda install path, so the container and a local conda env resolve to
# the same pinned versions. It still contains snakemake, which means this one
# image can run BOTH the legacy Snakemake pipeline and the Nextflow pipeline —
# that is what makes an output-level diff between them a fair comparison.
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

WORKDIR /data
