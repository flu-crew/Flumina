#### filter_flumut_subtype.R
####
#### Drop FluMut HA/NA markers when the reference subtype does not match the
#### numbering scheme used by those markers.
####
#### FluMut is an H5N1 tool, but its database is broader than that label suggests
#### and carries markers across many subtypes. Refusing to run it off-subtype
#### would discard most of its value. The internal genes are subtype-agnostic:
#### PB2, PB1, PA, NP, and NS markers are subtype-agnostic biology (PB2:K702R
#### raises polymerase activity in mammalian cells whatever the HA is).
####
#### HA and NA require special handling. The database protein names identify their
#### numbering schemes: HA1-5 uses H5 HA1 numbering, and NA-1 uses N1 numbering.
#### H3 and H5 HA1 differ in length and alignment, so
#### position 139 in one is not position 139 in the other. And the two have
#### diverged far enough that even a correctly mapped residue need not carry the
#### same meaning - a substitution that shifts receptor binding in H5 may do
#### nothing, or something else, in H3.
####
#### Filter on the protein prefix, never the Subtype column. Subtype records
#### which virus a finding was published in, not which numbering the position
#### uses: PB2:K702R is labelled H5N1 and is entirely valid on swine.
####
#### Determine subtype from segment names (e.g. A_HA_H3 / A_NA_N2, or A_HA_H5)
#### using two sources, in this order:
####
####   1. the reference FASTA's own segment names, and
####   2. the IRMA consensus contigs, whose headers IRMA writes with the subtype
####      it assigned (>A_HA_H5, >A_NA_N1), aggregated across samples.
####
#### The second source is necessary because a bare A_HA is unconfirmable and is
#### therefore treated as a mismatch. The bundled reference.fa is bare, but the
#### test_dataset samples are H5N1: IRMA assigns A_HA_H5 / A_NA_N1 to each, and
#### H5/N1 numbering is correct for those markers.
####
#### The reference takes precedence when it states a subtype, because that is an
#### explicit claim about the exact sequence FLUMUT_LOWFREQ screens in reference
#### coordinates. IRMA fills the gap only where the reference is silent - so this
#### is purely additive and cannot change a run that already resolved. Where the
#### two disagree nothing is silently reconciled: it is reported, because reads
#### assembling to a different subtype than the reference they were mapped to is
#### a problem with the run, not a labelling detail.
####
#### No BLAST or additional database is required: IRMA has already performed this
#### classification, and its result is available in a file staged by these
#### processes. The same assembly is trusted for all downstream analysis.
####
#### markers_all.tsv retains every row.
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

# IRMA's own call, one per sample, aggregated. A strict majority of samples
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

# Prefer the reference: it is an explicit claim about the exact sequence screened
# in reference coordinates. IRMA fills gaps but never overrides it.
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
    cat(sprintf("  %s: no majority across IRMA consensus (%s) - treated as unconfirmed\n",
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
# Record which source decided, so readers do not have to infer whether a subtype was
# stated by the reference or inferred from the assemblies.
cat(sprintf("  subtype source: HA from %s, NA from %s\n", ha.r$source, na.r$source),
    file = stderr())

# Publish the decision rather than leaving it only in the log.
#
# FluLens carries its own refSubtype() implementation of this rule against the
# reference names. Two copies can drift apart, so write the decision here and
# have FluLens read it instead of recomputing it.
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

# An absent markers.tsv ends the run only at this point, after the subtype file
# is written. The subtype is a property of the run, not of whether FluMut found
# anything. A check at the top would skip writing the file in exactly the cases
# where a reader most needs to know the screen was attempted.
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
  cat("  FLUMUT_KEEP_MISMATCHED_HA_NA is set - kept, but they are numbered for a\n",
      "  different subtype and should not be read as confirmed.\n", sep = "", file = stderr())
  quit(status = 0)
}
write.table(m[!drop, , drop = FALSE], file.path(outdir, "markers.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  dropped %d mismatched HA/NA row(s); %d remain (all rows kept in markers_all.tsv)\n",
            sum(drop), sum(!drop)), file = stderr())
