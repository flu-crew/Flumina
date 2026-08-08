#### filter_flumut_subtype.R
####
#### Drop FluMut HA/NA markers when the reference's subtype does not match the
#### one those markers are numbered for.
####
#### FluMut is an H5N1 tool. Its database is broader than that label suggests —
#### on the swine H3N2 run the 69 retained markers carried 16 subtype labels and
#### only 19 were H5N1 — so blanket-refusing to run it off-subtype throws away
#### most of its value. The internal genes are the value: PB2, PB1, PA, NP and
#### NS markers are subtype-agnostic biology (PB2:K702R raises polymerase
#### activity in mammalian cells whatever the HA is), and they were 50 of those
#### 69.
####
#### HA and NA are the problem, and the database's own protein names give it
#### away: HA1-5 is H5 HA1 numbering, NA-1 is N1 NA numbering. Two failures
#### compound off-subtype. H3 and H5 HA1 differ in length and alignment, so
#### position 139 in one is not position 139 in the other. And the two have
#### diverged far enough that even a correctly mapped residue need not carry the
#### same meaning — a substitution that shifts receptor binding in H5 may do
#### nothing, or something else, in H3.
####
#### Filter on the PROTEIN PREFIX, never the Subtype column. Subtype records
#### which virus a finding was published in, not which numbering the position
#### uses: PB2:K702R is labelled H5N1 and is entirely valid on swine.
####
#### Subtype comes from segment names — A_HA_H3 / A_NA_N2 on the swine reference,
#### A_HA_H5 on the cow data — from TWO sources, in this order:
####
####   1. the reference FASTA's own segment names, and
####   2. the IRMA consensus contigs, whose headers IRMA writes with the subtype
####      it assigned (>A_HA_H5, >A_NA_N1), aggregated across samples.
####
#### The second exists because a bare A_HA is unconfirmable, unconfirmed counts
#### as a mismatch, and the repo's OWN reference.fa is bare — so the conservative
#### branch was the common case rather than a corner one. It was also actively
#### wrong on the bundled test_dataset: those four samples are H5N1 (IRMA calls
#### every one A_HA_H5 / A_NA_N1) and their HA/NA markers were being dropped on
#### the one dataset here where H5/N1 numbering is exactly correct.
####
#### The reference still WINS when it states a subtype, because that is an
#### explicit claim about the exact sequence FLUMUT_LOWFREQ screens in reference
#### coordinates. IRMA fills the gap only where the reference is silent — so this
#### is purely additive and cannot change a run that already resolved. Where the
#### two disagree nothing is silently reconciled: it is reported, because reads
#### assembling to a different subtype than the reference they were mapped to is
#### a problem with the run, not a labelling detail.
####
#### No BLAST and no database: IRMA has already done this classification, its
#### answer is in a file these processes already stage, and we trust that same
#### assembly for everything else downstream.
####
#### Nothing is destroyed: markers_all.tsv already holds every row.
####
#### Usage:
####   Rscript filter_flumut_subtype.R <reference.fa> <markers.tsv> <keep> [outdir] [consensus_dir]
####     keep = TRUE keeps mismatched HA/NA rows (still annotated in the log)
####     consensus_dir = IRMA consensus contigs; optional, absent when IRMA is off

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript filter_flumut_subtype.R <reference.fa> <markers.tsv> <keep> [outdir]\n",
      file = stderr()); quit(status = 1)
}
ref.path <- args[1]; markers.path <- args[2]
keep <- toupper(args[3]) %in% c("TRUE", "T", "YES", "1")
outdir <- if (length(args) >= 4) args[4] else dirname(markers.path)
cons.dir <- if (length(args) >= 5) args[5] else NA_character_

seg_names_of <- function(path) {
  h <- grep("^>", readLines(path, warn = FALSE), value = TRUE)
  sub("^>", "", sub("[ \t].*$", "", h))
}

# "A_HA_H3" -> "H3"; "A_HA" -> NA (present but unconfirmable)
subtype_in <- function(seg.names, seg, letter) {
  hit <- grep(paste0("(^|_)", seg, "(_|$)"), seg.names, value = TRUE)
  if (!length(hit)) return(NA_character_)
  m <- regmatches(hit[1], regexpr(paste0(letter, "[0-9]+$"), hit[1]))
  if (!length(m)) NA_character_ else m
}

ref.names <- seg_names_of(ref.path)
ha.ref <- subtype_in(ref.names, "HA", "H")
na.ref <- subtype_in(ref.names, "NA", "N")

# IRMA's own call, one per sample, aggregated. A STRICT MAJORITY of the samples
# that produced a call is required, and the distribution is always printed: a
# cohort that genuinely splits across subtypes is a finding, not an error to be
# voted away, and silently taking the mode would hide it.
consensus_subtype <- function(dir, seg, letter) {
  out <- list(call = NA_character_, n = 0L, files = 0L, tab = NULL)
  if (is.na(dir) || !nzchar(dir) || !dir.exists(dir)) return(out)
  fa <- list.files(dir, pattern = "\\.(fa|fasta|fas)$", full.names = TRUE)
  out$files <- length(fa)
  if (!length(fa)) return(out)
  per <- vapply(fa, function(f)
    tryCatch(subtype_in(seg_names_of(f), seg, letter), error = function(e) NA_character_),
    character(1), USE.NAMES = FALSE)
  obs <- per[!is.na(per)]
  out$n <- length(obs)
  if (!out$n) return(out)
  tb <- sort(table(obs), decreasing = TRUE)
  out$tab <- tb
  if (tb[[1]] > out$n / 2) out$call <- names(tb)[1]
  out
}
ha.irma <- consensus_subtype(cons.dir, "HA", "H")
na.irma <- consensus_subtype(cons.dir, "NA", "N")

# Reference first: it is an explicit claim about the exact sequence screened in
# reference coordinates. IRMA only fills a gap, never overrides.
resolve <- function(ref, irma, gene) {
  if (!is.na(ref)) {
    if (!is.na(irma$call) && irma$call != ref)
      cat(sprintf("  WARNING: reference says %s=%s but IRMA assembled %s in %d of %d sample(s).\n",
                  gene, ref, irma$call, irma$n, irma$files),
          "           Using the reference. Reads assembling to a different subtype than the\n",
          "           reference they were mapped against is a problem with the run itself.\n",
          sep = "", file = stderr())
    return(list(value = ref, source = "reference name"))
  }
  if (!is.na(irma$call))
    return(list(value = irma$call,
                source = sprintf("IRMA consensus, %d of %d sample(s)", irma$n, irma$files)))
  if (!is.null(irma$tab) && irma$n > 0)
    cat(sprintf("  %s: no majority across IRMA consensus (%s) — treated as unconfirmed\n",
                gene, paste(names(irma$tab), irma$tab, sep = "=", collapse = " ")),
        file = stderr())
  list(value = NA_character_, source = "unconfirmed")
}
ha.r <- resolve(ha.ref, ha.irma, "HA"); na.r <- resolve(na.ref, na.irma, "NA")
ha <- ha.r$value; na <- na.r$value

ha.ok <- !is.na(ha) && ha == "H5"
na.ok <- !is.na(na) && na == "N1"
cat(sprintf("Reference subtype: HA=%s NA=%s -> HA markers %s, NA markers %s\n",
            if (is.na(ha)) "unconfirmed" else ha,
            if (is.na(na)) "unconfirmed" else na,
            if (ha.ok) "valid" else "MISMATCHED",
            if (na.ok) "valid" else "MISMATCHED"), file = stderr())
# Which source decided, so a reader never has to guess whether a subtype was
# stated by the reference or inferred from the assemblies.
cat(sprintf("  subtype source: HA from %s, NA from %s\n", ha.r$source, na.r$source),
    file = stderr())

# Publish the decision rather than leaving it in a log.
#
# FluLens carries its own refSubtype() applying "the same rule" against the
# reference names, and duplicated rules are exactly how this project has been
# bitten before — stats.multi disagreed with the panel it was meant to match for
# weeks because only one copy got fixed. The moment this script gained a second
# source, that copy became WRONG rather than merely duplicated: a bare-named
# H5N1 run now keeps its HA/NA markers here while FluLens would still call them
# unconfirmed. So the answer is written down and the app reads it instead of
# recomputing it.
#
# Tab-separated with a header, because everything else here is, and NA is
# written as the literal "unconfirmed" rather than an empty cell so it cannot be
# confused with a parse failure.
write.table(
  data.frame(gene   = c("HA", "NA"),
             subtype= c(if (is.na(ha)) "unconfirmed" else ha,
                        if (is.na(na)) "unconfirmed" else na),
             valid  = c(if (ha.ok) "TRUE" else "FALSE",
                        if (na.ok) "TRUE" else "FALSE"),
             source = c(ha.r$source, na.r$source),
             stringsAsFactors = FALSE),
  file.path(outdir, "subtype.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Only NOW does an absent markers.tsv end the run. The subtype is a property of
# the RUN, not of whether FluMut happened to find anything, and this check used
# to sit at the top — which would have skipped writing the file in exactly the
# cases where a reader most needs to know the screen was attempted.
if (!file.exists(markers.path) || file.info(markers.path)$size == 0) {
  cat("No markers to filter.\n", file = stderr()); quit(status = 0)
}

m <- read.delim(markers.path, sep = "\t", check.names = FALSE,
                stringsAsFactors = FALSE, quote = "")
if (!"Marker" %in% names(m)) {
  cat("markers.tsv has no Marker column; nothing to do.\n", file = stderr()); quit(status = 0)
}
protein <- sub(":.*$", "", m$Marker)
is.ha <- grepl("^HA", protein); is.na.seg <- grepl("^NA", protein)
drop <- (is.ha & !ha.ok) | (is.na.seg & !na.ok)

cat(sprintf("  HA marker rows: %d, NA marker rows: %d, mismatched: %d of %d\n",
            sum(is.ha), sum(is.na.seg), sum(drop), nrow(m)), file = stderr())

if (!any(drop)) {
  cat("  nothing to drop\n", file = stderr()); quit(status = 0)
}
if (keep) {
  cat("  FLUMUT_KEEP_MISMATCHED_HA_NA is set — kept, but they are numbered for a\n",
      "  different subtype and should not be read as confirmed.\n", sep = "", file = stderr())
  quit(status = 0)
}
write.table(m[!drop, , drop = FALSE], file.path(outdir, "markers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  dropped %d mismatched HA/NA row(s); %d remain (all rows kept in markers_all.tsv)\n",
            sum(drop), sum(!drop)), file = stderr())
