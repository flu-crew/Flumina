# IRMA parameter file — example
#
# Flumina passes this to IRMA with --external-config. A file named
# irma_config.sh in your working directory is picked up automatically; point at
# a different one with `flumina -x /path/to/irma_config.sh`.
#
# You very likely do not need this file at all. IRMA's defaults are sensible,
# and IRMA 1.3.5 works out how many cores it may use on its own — Flumina also
# tells it exactly what the scheduler granted. Everything below is commented
# out on purpose; uncomment only what you actually mean to change.


# Temporary directory IRMA writes its intermediates to.
#
# WARNING: IRMA does NOT create this directory. It builds its working path
# directly from the value ("$TMP"/user/IRMAvX/run-token), so a path that does
# not exist makes the match stage silently produce nothing: reads are counted,
# then every table and consensus comes out empty and IRMA still exits 0.
#
# Leave this unset unless you have a reason. If you do set it, use an absolute
# path that already exists, such as fast local scratch on a compute node.
#TMP=/scratch/myuser/irma_tmp


# How many processes IRMA may run at once.
#
# Leave unset. IRMA 1.3.5 detects the cores available to it, and Flumina exports
# LOCAL_PROCS_OVERRIDE from the CPUs the scheduler actually granted the job, so
# these are already correct. Setting them by hand is how you end up
# oversubscribing a node or leaving most of an allocation idle.
#SINGLE_LOCAL_PROC=2
#DOUBLE_LOCAL_PROC=1
