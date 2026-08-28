#!/bin/bash
#
# Example cluster job script for running Flumina — command-line arguments version.
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
#     sbatch job_script_example_arguments.sh   (or qsub on PBS)
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

# ---- On PBS/Torque or PBS Pro instead? -------------------------------------
# Replace the #SBATCH block above with this, set PROFILE below to pbs (or
# pbspro), and submit with qsub rather than sbatch:
#
#     #PBS -N Flumina
#     #PBS -q priority
#     #PBS -A nadc_iav
#     #PBS -l select=1:ncpus=2:mem=8gb
#     #PBS -l walltime=168:00:00
#     #PBS -m abe
#     #PBS -M your.email@example.com
#
# Two differences that catch people out: PBS starts the job in your HOME rather
# than where you submitted from, and the job id is $PBS_JOBID, not
# $SLURM_JOB_ID. Both are handled just below. The flumina command is identical.
# ----------------------------------------------------------------------------

date

# Works under either scheduler: PBS needs an explicit cd, Slurm does not.
cd "${PBS_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$PWD}}" || exit 1

# Job id under whichever scheduler is running this, used to keep concurrent
# runs from sharing a work directory.
JOB_ID="${SLURM_JOB_ID:-${PBS_JOBID:-$$}}"

#############################################
#### Settings
#############################################

# Where your paired fastq.gz reads live
READS="raw_reads"

# Where results should be written
OUTPUT="Flumina_results"

# Nextflow's intermediate files. Only set this if the directory you submit from
# is a poor place for it — home or project space, with quotas and slow shared
# storage. It is written to constantly and grows to many times the size of the
# final results, and with per-step job submission every compute node must be
# able to see it. If you already submit from fast scratch, delete these two
# lines and let it default to ./work alongside your results.
#
# Deliberately NOT scoped to the job id: the work directory IS the resume cache,
# so a per-job path would silently make -R useless — every resubmission would
# start from an empty directory and redo everything. Naming it after the output
# keeps resume working across submissions while still separating unrelated runs.
WORKDIR="/scratch/${USER}/flumina_work/$(basename "$OUTPUT")"

# Which scheduler to submit the per-step jobs to, paired with apptainer so each
# of them runs inside the container. Options: slurm, pbs (Torque/OpenPBS),
# pbspro, sge, lsf.
PROFILE="slurm,apptainer"

# Scheduler settings for the jobs Flumina submits. These are the flu-crew
# defaults and WILL be rejected on any other cluster — change them. The account
# string is passed through untouched, so use your scheduler's own syntax.
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
FLUMINA_DIR="${PWD}"
if ! command -v flumina &> /dev/null; then
    export PATH="${FLUMINA_DIR}/Scripts:${PATH}"
fi

flumina --version


#############################################
#### Run Flumina
#############################################

# The PROFILE set above is what makes this parallel: the scheduler half submits
# each step as its own job, `apptainer` runs each of those inside the container.
# Nextflow converts the image once into a shared cache before launching
# anything, so the submitted jobs do not each race to build their own copy.

mkdir -p "$WORKDIR"

flumina \
    -i "$READS" \
    -o "$OUTPUT" \
    -p "$PROFILE" \
    -w "$WORKDIR" \
    -Q "$QUEUE" \
    -A "$ACCOUNT" \
    -j "$MAX_JOBS" \
    -t "$THREADS_PER_STEP" \
    -M "$MEMORY_PER_STEP"

# Nextflow keeps every intermediate file in the work directory, so it is worth
# reclaiming once you are done — but NOT automatically. Deleting it here throws
# away the resume cache, so a run that fails at hour twenty has to start over
# from nothing. Clear it by hand when you are satisfied with the results:
#
#     rm -rf "$WORKDIR"

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
#     READS="/project/sequencing/run47/fastq"
#     OUTPUT="Flumina_results"
#     WORKDIR="/scratch/${USER}/flumina_${JOB_ID}"
#
#     # Apptainer mounts your home and current directory automatically but
#     # nothing else, so reads or scratch anywhere else are invisible inside the
#     # container. Rather than work out --bind flags by hand, take the top-level
#     # directory of each path and bind those. APPTAINER_BIND is picked up by
#     # `apptainer run` without changing the command line.
#     binds=""
#     for p in "$READS" "$OUTPUT" "$WORKDIR"; do
#         case "$p" in
#           /*) top="/$(printf '%s' "${p#/}" | cut -d/ -f1)"
#               case ",$binds," in *",$top,"*) ;; *) binds="${binds:+$binds,}$top" ;; esac ;;
#         esac
#     done
#     [ -n "$binds" ] && export APPTAINER_BIND="$binds"
#
#     mkdir -p "$WORKDIR"
#     apptainer run flumina_latest.sif \
#         -i "$READS" -o "$OUTPUT" \
#         -t "${SLURM_CPUS_PER_TASK:-${NCPUS:-24}}" -M 100.GB \
#         -w "$WORKDIR"
#
# Steps then run in parallel only up to the cores of that one node, so this is
# considerably slower for large sample counts — but it needs nothing installed
# beyond Apptainer.
