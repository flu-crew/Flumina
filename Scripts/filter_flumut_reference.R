#### filter_flumut_reference.R
####
#### Remove the reference's own FluMut findings from a set of sample results.
####
#### FluMut reports every marker carried by a sequence, including markers already
#### present in the reference. Such markers occur in every sample by construction
#### and provide no sample-specific information.
####
#### Usage:
####   Rscript filter_flumut_reference.R <ref_markers> <ref_mutations> \
####       <markers> <mutations> <literature> <outdir>
####
#### Write the following files to <outdir>:
####   markers.tsv        sample rows not present in the reference
####   mutations.tsv      columns where at least one sample differs from the reference
####   literature.tsv     literature for the retained markers
####   *_all.tsv          the unfiltered originals, always kept
####   reference_*.tsv    markers carried by the reference, so each removal is
####                      traceable
####
#### Preserve every input under *_all.tsv.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat("Usage: Rscript filter_flumut_reference.R <ref_markers> <ref_mutations>",
      "<markers> <mutations> <literature> <outdir>\n", file = stderr())
  quit(status = 1)
}
ref.markers.p   <- args[1]
ref.mutations.p <- args[2]
markers.p       <- args[3]
mutations.p     <- args[4]
literature.p    <- args[5]
outdir          <- args[6]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# `colClasses = "character"` is required for correctness. These tables contain
# residues and names, but read.table still infers a type for each column. Because
# the reference tables contain one row, a column containing only "T" or "F" can
# be inferred as logical and written as TRUE/FALSE.
#
# This also breaks the invariance test below: `rv` becomes "TRUE" while samples
# contain "T", so `all(sv == rv)` is never true and an invariant column is kept as
# informative.
rd <- function(p) {
  if (!file.exists(p) || file.info(p)$size == 0) return(NULL)
  read.table(p, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, na.strings = "",
             colClasses = "character")
}
wr <- function(d, p) write.table(d, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

ref.markers   <- rd(ref.markers.p)
ref.mutations <- rd(ref.mutations.p)
markers       <- rd(markers.p)
mutations     <- rd(mutations.p)
literature    <- rd(literature.p)

# Preserve the original files before filtering.
if (!is.null(markers))    wr(markers,    file.path(outdir, "markers_all.tsv"))
if (!is.null(mutations))  wr(mutations,  file.path(outdir, "mutations_all.tsv"))
if (!is.null(literature)) wr(literature, file.path(outdir, "literature_all.tsv"))
if (!is.null(ref.markers))   wr(ref.markers,   file.path(outdir, "reference_markers.tsv"))
if (!is.null(ref.mutations)) wr(ref.mutations, file.path(outdir, "reference_mutations.tsv"))

#############################################
#### markers.tsv — remove rows also present in the reference
#############################################
# Long format: Sample | Marker | Mutations in your sample | Effect | Subtype | Literature
# The key is the complete row except Sample and Literature. "Mutations in your
# sample" is included because under --relaxed a marker can be reported from
# a partial match, and a sample matching more mutations than the reference did
# is a different finding even though the Marker string is the same.
if (!is.null(markers) && nrow(markers) > 0) {
  key.cols <- intersect(c("Marker", "Mutations in your sample", "Effect", "Subtype"),
                        colnames(markers))
  mk_key <- function(d) do.call(paste, c(d[, key.cols, drop = FALSE], sep = "\r"))

  if (!is.null(ref.markers) && nrow(ref.markers) > 0 &&
      all(key.cols %in% colnames(ref.markers))) {
    ref.keys <- unique(mk_key(ref.markers))
    keep     <- !(mk_key(markers) %in% ref.keys)
  } else {
    ref.keys <- character(0)
    keep     <- rep(TRUE, nrow(markers))
    message("No reference markers to subtract; markers.tsv left unfiltered.")
  }

  out.markers <- markers[keep, , drop = FALSE]
  wr(out.markers, file.path(outdir, "markers.tsv"))
  cat(sprintf("markers.tsv:   %d rows -> %d (%d shared with the reference removed)\n",
              nrow(markers), nrow(out.markers), sum(!keep)))

  # literature.tsv follows the markers that survived filtering.
  if (!is.null(literature) && nrow(literature) > 0) {
    lit.key <- intersect(c("Marker", "Effect", "Subtype"), colnames(literature))
    if (length(lit.key) > 0 && all(lit.key %in% colnames(out.markers))) {
      kept.lit <- unique(do.call(paste, c(out.markers[, lit.key, drop = FALSE], sep = "\r")))
      lit.keep <- do.call(paste, c(literature[, lit.key, drop = FALSE], sep = "\r")) %in% kept.lit
      wr(literature[lit.keep, , drop = FALSE], file.path(outdir, "literature.tsv"))
      cat(sprintf("literature.tsv: %d rows -> %d\n", nrow(literature), sum(lit.keep)))
    } else {
      wr(literature, file.path(outdir, "literature.tsv"))
    }
  } else if (!is.null(literature)) {
    wr(literature, file.path(outdir, "literature.tsv"))
  }
}

#############################################
#### mutations.tsv — remove columns with no sample-level variation
#############################################
# Wide format: Sample | <one column per mutation>, cell = the residue found.
#
# A column is dropped only when every sample carries the reference residue.
# Dropping by "the reference has this marker" would discard reversions: a sample
# that loses a reference marker produces no marker row, so markers.tsv cannot
# show it. The wide table is therefore the only place that signal exists.
if (!is.null(mutations) && nrow(mutations) > 0) {
  sample.col <- colnames(mutations)[1]
  mut.cols   <- setdiff(colnames(mutations), sample.col)

  if (!is.null(ref.mutations) && nrow(ref.mutations) > 0) {
    ref.row <- ref.mutations[1, , drop = FALSE]
    informative <- vapply(mut.cols, function(cl) {
      if (!cl %in% colnames(ref.row)) return(TRUE)   # reference never saw it: keep
      rv <- as.character(ref.row[[cl]])
      sv <- unique(as.character(mutations[[cl]]))
      # Drop missing values before comparing. A sample with no data for this
      # marker is not evidence of a difference. An NA left in makes `sv == rv`
      # return NA, which `all()` propagates and which then indexes the column
      # list as NA ("undefined columns selected"). Samples can miss whole
      # segments, so this is common.
      sv <- sv[!is.na(sv)]
      if (length(sv) == 0L) return(FALSE)            # nothing observed: no signal
      if (is.na(rv)) return(TRUE)                    # reference blank, samples not
      !all(sv == rv)
    }, logical(1))
  } else {
    informative <- rep(TRUE, length(mut.cols))
    message("No reference mutations to compare; mutations.tsv left unfiltered.")
  }

  informative[is.na(informative)] <- TRUE           # never index columns by NA
  out.mutations <- mutations[, c(sample.col, mut.cols[informative]), drop = FALSE]
  wr(out.mutations, file.path(outdir, "mutations.tsv"))
  cat(sprintf("mutations.tsv: %d mutation columns -> %d (%d invariant vs the reference removed)\n",
              length(mut.cols), sum(informative), sum(!informative)))
}

cat("Reference findings kept in reference_markers.tsv / reference_mutations.tsv;",
    "unfiltered output in *_all.tsv\n")
