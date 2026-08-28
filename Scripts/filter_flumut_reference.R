#### filter_flumut_reference.R
####
#### Remove the reference's own FluMut findings from a set of sample results.
####
#### FluMut reports every marker a sequence carries, including the ones the
#### reference already carries. Those are present in every sample by
#### construction, so they say nothing about any sample. On the swine WGS run
#### the reference alone accounts for 84 marker rows, and across 30 samples
#### 2,514 of 2,515 reported rows were identical to it — one row in 2,515
#### carried information.
####
#### Usage:
####   Rscript filter_flumut_reference.R <ref_markers> <ref_mutations> \
####       <markers> <mutations> <literature> <outdir>
####
#### Writes into <outdir>:
####   markers.tsv        sample rows the reference does NOT have
####   mutations.tsv      columns where at least one sample differs from the reference
####   literature.tsv     literature for the retained markers
####   *_all.tsv          the unfiltered originals, always kept
####   reference_*.tsv    what the reference itself carries, so a removal can
####                      always be explained
####
#### Nothing is destroyed: every input is preserved under *_all.tsv.

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

# colClasses = "character" is load-bearing, not tidiness. Every cell in these
# tables is a residue or a name, but read.table still type-GUESSES per column,
# and the reference tables have exactly one row -- so a column whose only value
# is "T" (threonine) or "F" (phenylalanine) is guessed as logical and written
# back out as TRUE/FALSE. On the swine WGS run that corrupted 8 of 59 residues
# in reference_mutations.tsv.
#
# It also silently breaks the invariance test below: rv becomes "TRUE" while the
# samples say "T", so `all(sv == rv)` is never true and the column is kept as
# informative even when every sample matches the reference. That has not changed
# a result yet -- on this run all 8 columns vary on their own merits -- but it is
# a live trap the moment a T/F column is genuinely invariant.
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

# Preserve the originals first, unconditionally.
if (!is.null(markers))    wr(markers,    file.path(outdir, "markers_all.tsv"))
if (!is.null(mutations))  wr(mutations,  file.path(outdir, "mutations_all.tsv"))
if (!is.null(literature)) wr(literature, file.path(outdir, "literature_all.tsv"))
if (!is.null(ref.markers))   wr(ref.markers,   file.path(outdir, "reference_markers.tsv"))
if (!is.null(ref.mutations)) wr(ref.mutations, file.path(outdir, "reference_mutations.tsv"))

#############################################
#### markers.tsv — drop rows the reference also has
#############################################
# Long format: Sample | Marker | Mutations in your sample | Effect | Subtype | Literature
# The key is the whole row EXCEPT Sample and Literature. "Mutations in your
# sample" is included on purpose: under --relaxed a marker can be reported from
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

  # literature.tsv follows whichever markers survived
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
#### mutations.tsv — drop columns no sample varies at
#############################################
# Wide format: Sample | <one column per mutation>, cell = the residue found.
#
# A column is dropped only when EVERY sample carries the reference residue.
# Dropping by "the reference has this marker" instead would discard reversions,
# which is the one thing markers.tsv can never show — a sample that LOSES a
# reference marker simply produces no marker row, so the wide table is the only
# place that signal exists. On the swine WGS run this rule cut 59 columns to 4,
# and one of the four was exactly such a reversion (NA-1:S364N, reference N,
# one sample Q).
if (!is.null(mutations) && nrow(mutations) > 0) {
  sample.col <- colnames(mutations)[1]
  mut.cols   <- setdiff(colnames(mutations), sample.col)

  if (!is.null(ref.mutations) && nrow(ref.mutations) > 0) {
    ref.row <- ref.mutations[1, , drop = FALSE]
    informative <- vapply(mut.cols, function(cl) {
      if (!cl %in% colnames(ref.row)) return(TRUE)   # reference never saw it: keep
      rv <- as.character(ref.row[[cl]])
      sv <- unique(as.character(mutations[[cl]]))
      # Drop missing values BEFORE comparing. A sample with no data for this
      # marker is not evidence of a difference, and leaving the NA in makes
      # `sv == rv` return NA, which `all()` propagates and which then indexes
      # the column list as NA ("undefined columns selected"). The swine WGS run
      # has samples missing whole segments, so this is the normal case, not an
      # edge one.
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
