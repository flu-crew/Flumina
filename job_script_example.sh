#!/bin/bash
#
# Example SLURM job script for running Flumina across a cluster.
#
# This is the FAST mode. Nextflow runs here, in this small job, and submits every
# pipeline step as its own SLURM job — so independent samples and independent
# steps run at the same time instead of one after another. With 100 samples the
# per-sample chain runs ~100-wide, which is where Flumina's speed comes from.
#
# Each of those submitted jobs runs inside the Flumina Apptainer container, so
# you still get exactly one pinned software environment. Nextflow itself has to
# run on the host rather than in the container, because it submits work with
# sbatch and that binary is not inside the image.
#
# This job is deliberately TINY — 2 cores is plenty. It spends its life waiting
# on the jobs it submits, so do not request a big node for it, and do give it a
# long wall time, since it must outlive everything it launches.
#
# Submit with:
#     sbatch job_script_example.sh
#
# EDIT BEFORE SUBMITTING: the account/partition lines, the email address, the
# QUEUE/ACCOUNT settings, and the paths under "Settings".

#SBATCH --job-name=Flumina
#SBATCH -p priority --qos=vpru
#SBATCH -A nadc_iav
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 168:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@example.com
#SBATCH -o "flumina.%j.out"
#SBATCH -e "flumina.%j.err"

date

#############################################
#### Settings
#############################################

# Where your paired fastq.gz reads live
READS="raw_reads"

# Where results should be written
OUTPUT="Flumina_results"

# Nextflow's intermediate files. This MUST be on a filesystem every compute node
# can see, because the submitted jobs all read and write here. Keep it on
# scratch, not home: it is written to constantly and grows large.
WORKDIR="/scratch/${USER}/flumina_${SLURM_JOB_ID}"

# Scheduler settings for the jobs Flumina submits. These are the flu-crew
# defaults and WILL be rejected on any other cluster — change them.
QUEUE="priority"
ACCOUNT="--qos=vpru -A nadc_iav"

# How many jobs to keep in flight at once. This is the main throughput dial:
# raise it to go wider if your allocation allows, lower it to be a better
# neighbour on a busy shared queue.
MAX_JOBS=100

# Resources for each individual step. These apply per job, not to the whole run.
THREADS_PER_STEP=8
MEMORY_PER_STEP="32.GB"


#############################################
#### Load Nextflow and Apptainer
#############################################

# Nextflow runs on the host here, so it must be available as a module or on your
# PATH. Most clusters provide both of these; adjust the names to match yours. If
# Apptainer is still called Singularity, load that instead — Nextflow handles it.
module load nextflow
module load apptainer

# Add Flumina's Scripts folder to PATH if `flumina` is not already there. Point
# FLUMINA_DIR at wherever you cloned this repository.
FLUMINA_DIR="${SLURM_SUBMIT_DIR}"
if ! command -v flumina &> /dev/null; then
    export PATH="${FLUMINA_DIR}/Scripts:${PATH}"
fi

flumina --version


#############################################
#### Run Flumina
#############################################

# -p slurm,apptainer is what makes this parallel: `slurm` submits each step as
# its own job, `apptainer` runs each of those inside the container. Nextflow
# converts the image once into a shared cache before launching anything, so the
# submitted jobs do not each race to build their own copy.

mkdir -p "$WORKDIR"

flumina \
    -i "$READS" \
    -o "$OUTPUT" \
    -p slurm,apptainer \
    -w "$WORKDIR" \
    -Q "$QUEUE" \
    -A "$ACCOUNT" \
    -j "$MAX_JOBS" \
    -t "$THREADS_PER_STEP" \
    -M "$MEMORY_PER_STEP"

# Nextflow keeps every intermediate file in the work directory. Once results are
# published it is only needed for -R (resume), so clear it to release the space.
# Comment this out while you are still iterating on a run.
rm -rf "$WORKDIR"

date

# End of file
#
#############################################
#### The simpler alternative
#############################################
#
# If your cluster has no Nextflow module and you would rather not install it, or
# you are only running a handful of samples, everything can instead run inside a
# single allocation with no host-side Nextflow at all:
#
#     #SBATCH --cpus-per-task=24
#     #SBATCH --mem=100G
#
#     module load apptainer
#     [ -f flumina_latest.sif ] || apptainer pull flumina_latest.sif docker://chutter/flumina
#
#     apptainer run flumina_latest.sif \
#         -i raw_reads -o Flumina_results \
#         -t "${SLURM_CPUS_PER_TASK}" -M 100.GB \
#         -w "/scratch/${USER}/flumina_${SLURM_JOB_ID}"
#
# Steps then run in parallel only up to the cores of that one node, so this is
# considerably slower for large sample counts — but it needs nothing installed
# beyond Apptainer.
