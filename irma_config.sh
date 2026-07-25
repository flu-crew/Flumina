# IRMA parameter file — example
#
# IRMA reads this from its working directory on its own. Flumina picks up a file
# named irma_config.sh in your working directory automatically, or point at a
# different one with `flumina -x /path/to/irma_config.sh`.
#
# Any standard IRMA parameter can be set here; the three below are the ones most
# worth setting. See the IRMA documentation for the rest.

# Temporary directory IRMA writes to. IRMA's own default can land somewhere slow
# or unwritable, so set it explicitly. This relative path resolves inside the
# per-task working directory and is cleaned up with the rest of the run; give an
# absolute path instead if your cluster provides a fast local scratch disk.
TMP=./irma_tmp

# Should match the CPUs available to a single IRMA job (see `flumina -t`).
SINGLE_LOCAL_PROC=2
DOUBLE_LOCAL_PROC=1
