# Example IRMA parameter file.
#
# Flumina passes this to IRMA with --external-config. A file named
# irma_config.sh in your working directory is picked up automatically; point at
# a different one with `flumina -x /path/to/irma_config.sh`.
#
# This file is usually unnecessary. IRMA's defaults are appropriate, and IRMA
# 1.3.5 determines its own core usage; Flumina also tells it exactly what the
# scheduler granted. The settings below are commented out intentionally; enable
# only those that must be changed.


# Temporary directory for IRMA intermediate files.
#
# Warning: IRMA does not create this directory. It builds its working path
# directly from the value ("$TMP"/user/IRMAvX/run-token), so a path that does
# does not exist, the match stage produces no output: reads are counted, every
# table and consensus is empty, and IRMA still exits 0.
#
# Leave this unset unless you have a reason. If you do set it, use an absolute
# path that already exists, such as fast local scratch on a compute node.
#TMP=/scratch/myuser/irma_tmp


# Maximum number of concurrent IRMA processes.
#
# Leave unset. IRMA 1.3.5 detects available cores, and Flumina exports
# LOCAL_PROCS_OVERRIDE from the CPUs the scheduler actually granted the job, so
# the scheduler allocation. Setting them manually can oversubscribe a node or
# leave part of the allocation idle.
#SINGLE_LOCAL_PROC=2
#DOUBLE_LOCAL_PROC=1

# Minimum read length to pass IRMA's pre-processing quality control.
#
# IRMA's default for the FLU module is 125bp. If your sequencing run used 
# shorter reads (e.g., 100bp), IRMA will discard them all unless you lower 
# this threshold.
MIN_LEN=80
