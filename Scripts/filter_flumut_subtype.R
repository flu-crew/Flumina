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
#### Subtype is read from the reference segment names — A_HA_H3 / A_NA_N2 on the
#### swine reference, A_HA_H5 on the cow data. A bare A_HA carries no subtype, so
#### it CANNOT be confirmed, and unconfirmed counts as a mismatch rather than a
#### pass. The repo's own reference.fa is bare, so this is the common case, not
#### a corner one.
####
#### Nothing is destroyed: markers_all.tsv already holds every row.
####
#### Usage:
####   Rscript filter_flumut_subtype.R <reference.fa> <markers.tsv> <keep> [outdir]
####     keep = TRUE keeps mismatched HA/NA rows (still annotated in the log)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript filter_flumut_subtype.R <reference.fa> <markers.tsv> <keep> [outdir]\n",
      file = stderr()); quit(status = 1)
}
ref.path <- args[1]; markers.path <- args[2]
keep <- toupper(args[3]) %in% c("TRUE", "T", "YES", "1")
outdir <- if (length(args) >= 4) args[4] else dirname(markers.path)

if (!file.exists(markers.path) || file.info(markers.path)$size == 0) {
  cat("No markers to filter.\n", file = stderr()); quit(status = 0)
}

headers <- grep("^>", readLines(ref.path, warn = FALSE), value = TRUE)
seg.names <- sub("^>", "", sub("[ \t].*$", "", headers))

# "A_HA_H3" -> "H3"; "A_HA" -> NA (present but unconfirmable)
subtype_of <- function(seg, letter) {
  hit <- grep(paste0("(^|_)", seg, "(_|$)"), seg.names, value = TRUE)
  if (!length(hit)) return(NA_character_)
  m <- regmatches(hit[1], regexpr(paste0(letter, "[0-9]+$"), hit[1]))
  if (!length(m)) NA_character_ else m
}
ha <- subtype_of("HA", "H"); na <- subtype_of("NA", "N")

ha.ok <- !is.na(ha) && ha == "H5"
na.ok <- !is.na(na) && na == "N1"
cat(sprintf("Reference subtype: HA=%s NA=%s -> HA markers %s, NA markers %s\n",
            if (is.na(ha)) "unconfirmed" else ha,
            if (is.na(na)) "unconfirmed" else na,
            if (ha.ok) "valid" else "MISMATCHED",
            if (na.ok) "valid" else "MISMATCHED"), file = stderr())

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
