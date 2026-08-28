#!/bin/bash
#
# Example SLURM job script — configuration-file version.
#
# Everything about the analysis lives in config.cfg, so this script only has to
# load the software and start the run. Use this when you have more than a
# handful of non-default settings, or when you want a single file that records
# exactly how a run was configured. For the version that sets everything as
# command-line arguments instead, see job_script_example_arguments.sh.
#
# Setup, once:
#     module load apptainer
#     apptainer pull -F docker://chutter/flumina
#     apptainer run flumina_latest.sif --export ./flumina
#     cp flumina/config.cfg .          # then edit it for your run
#
# Submit with:
#     sbatch job_script_example_config.sh
#
# EDIT BEFORE SUBMITTING: the account/partition and email lines below, and
# config.cfg itself.

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
# Replace the #SBATCH block above with this, set PROFILE="pbs,apptainer" (or
# pbspro) in config.cfg, and submit with qsub rather than sbatch:
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
# than where you submitted from, so add `cd "$PBS_O_WORKDIR"` below; and the job
# id is $PBS_JOBID, not $SLURM_JOB_ID. Both are handled automatically further
# down. Everything else — the flumina command itself — is identical.
# ----------------------------------------------------------------------------

date

# Works under either scheduler: PBS needs an explicit cd, Slurm does not.
cd "${PBS_O_WORKDIR:-${SLURM_SUBMIT_DIR:-$PWD}}" || exit 1

# Job id under whichever scheduler is running this, used to keep concurrent
# runs from sharing a work directory.
JOB_ID="${SLURM_JOB_ID:-${PBS_JOBID:-$$}}"

# This job only runs Nextflow, which then submits each pipeline step as its own
# job and waits. Two cores is plenty; the long wall time above matters, because
# it has to outlive everything it launches.
module load nextflow apptainer

# Where the pipeline was exported to (see the setup block above)
export PATH="${PWD}/flumina/Scripts:${PATH}"

flumina --version

# Every setting comes from config.cfg, including the input and output paths,
# the profile, and the cluster queue. To submit steps as separate jobs, set
# PROFILE="slurm,apptainer" in there.
#
flumina -c config.cfg

# Nothing needs overriding on the command line here: the work directory defaults
# to ./work beside your results, which is right whenever you submit from fast
# scratch. Set WORK_DIRECTORY in config.cfg only if you submit from home or
# project space, where a constantly-written, fast-growing directory does not
# belong.
#
# Do not scope that path to the job id: the work directory IS the resume cache,
# so a per-job path would quietly make -R useless, and every resubmission would
# redo everything from scratch.
#
# Reclaim the space by hand once you are satisfied with the results, rather than
# automatically here — deleting it is what makes a failed run unresumable:
#
#     rm -rf work

date

# End of file
