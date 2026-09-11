#!/usr/bin/env Rscript

#### hana_transferable.R
####
#### Identify which off-subtype HA/NA markers are untransferable and which are
#### only labelled that way.
####
#### FLUMUT_KEEP_MISMATCHED_HA_NA drops every HA/NA marker when the reference's
#### subtype does not match FluMut's H5/N1 numbering. The case for that is real
#### but it is an average: FluMut's reference and this run's carry different
#### wild-type residues at most HA1 positions, so a marker phrased as a deviation
#### from the H5 residue is read against a protein that never had it.
####
#### The answer is known for each marker. flumut_position_map.tsv records, for
#### every position it places, whether the two references carry the same residue.
#### A marker on an identical residue is a deviation from the same starting point,
#### so its H5 provenance is a label rather than a confound.
####
#### Read a finished run and write one row per dropped HA/NA marker with that
#### verdict.
####
#### Rows are counted AFTER reference subtraction -- markers_all.tsv minus
#### reference_markers.tsv. markers_all.tsv alone inflates the answer with markers
#### this run's own reference carries, which appear in nearly every sample.
####
#### usage: Rscript hana_transferable.R <run_dir> [out.tsv]

MARK = "^([A-Za-z0-9-]+):([A-Z])([0-9]+)([A-Z])$"

# Return one column (1-based) of a headed TSV, skipping the header, or an empty
# vector when the file is absent.
rows.of = function(path, col = 2L) {
  if (!file.exists(path)) return(character(0))
  lines = readLines(path)
  if (length(lines) <= 1L) return(character(0))
  vals = character(0)
  for (line in lines[-1]) {
    f = strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) >= col) vals = c(vals, f[col])
  }
  vals
}

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  cat("Usage: Rscript hana_transferable.R <run_dir> [out.tsv]\n", file = stderr())
  quit(status = 1)
}
run = args[1]
fm = file.path(run, "variant_analysis", "flumut")
out = if (length(args) >= 2L) args[2] else file.path(fm, "hana_transferable.tsv")

# Position map: label -> flumut_pos -> (flumut_aa, ref_aa, identical).
pmap = file.path(fm, "flumut_position_map.tsv")
if (!file.exists(pmap)) {
  cat("no position map at ", pmap, " -- run FLUMUT_POSITION_MAP first\n",
      sep = "", file = stderr())
  quit(status = 1)
}
pos = list()
plines = readLines(pmap)
for (line in plines[-1]) {
  f = strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(f) < 8L) next
  label = f[1]
  if (is.null(pos[[label]])) pos[[label]] = list()
  pos[[label]][[f[3]]] = c(f[4], f[7], f[8])   # flumut_aa, ref_aa, identical
}

refm = unique(rows.of(file.path(fm, "reference_markers.tsv")))
kept = unique(rows.of(file.path(fm, "markers.tsv")))

# Count markers in markers_all.tsv that the reference does not carry, keeping
# first-appearance order for ties.
all.markers = rows.of(file.path(fm, "markers_all.tsv"))
all.markers = all.markers[!(all.markers %in% refm)]
uniq = unique(all.markers)
count.vec = as.integer(table(factor(all.markers, levels = uniq)))
ord = order(-count.vec)
markers.sorted = uniq[ord]
counts.sorted = count.vec[ord]

recs = list()
same = 0L
for (m in seq_along(markers.sorted)) {
  marker = markers.sorted[m]
  n = counts.sorted[m]
  for (tok in trimws(strsplit(marker, ",", fixed = TRUE)[[1]])) {
    g = regmatches(tok, regexec(MARK, tok))[[1]]
    if (length(g) == 0L) next
    label = g[2]; wt = g[3]; p = g[4]
    if (!substr(label, 1L, 2L) %in% c("HA", "NA")) next
    is.kept = if (marker %in% kept) "yes" else "no"
    rec = if (!is.null(pos[[label]])) pos[[label]][[p]] else NULL
    if (is.null(rec)) {
      recs[[length(recs) + 1L]] = list(marker, tok, label, p, "", "",
                                       "unplaced", n, is.kept, "")
      next
    }
    fm.aa = rec[1]; ref.aa = rec[2]; ident = rec[3]
    if (ident == "yes") same = same + 1L
    # The marker's wild-type should match FluMut's residue. When it does not, the
    # marker is phrased against a different reference residue, so the verdict is
    # less reliable.
    note = if (wt == fm.aa) "" else paste0("marker wild-type ", wt)
    recs[[length(recs) + 1L]] = list(marker, tok, label, p, fm.aa, ref.aa,
                                     ident, n, is.kept, note)
  }
}

header = paste("marker", "component", "label", "flumut_pos", "flumut_aa",
               "ref_aa", "wt_identical", "sample_rows", "kept", "note",
               sep = "\t")
body = vapply(recs, function(r) paste(unlist(r), collapse = "\t"), character(1))
writeLines(c(header, body), out)

# One row per component, so a compound marker appears twice. Sum sample rows over
# distinct markers to avoid double-counting compound markers.
distinct.markers = unique(vapply(recs, function(r) r[[1]], character(1)))
distinct.sum = sum(counts.sorted[match(distinct.markers, markers.sorted)])
cat("HA/NA after reference subtraction: ", length(distinct.markers),
    " distinct markers, ", distinct.sum, " sample-rows, ", length(recs),
    " components\n", sep = "")
cat("  components on an IDENTICAL wild-type residue (transferable): ", same, "\n",
    sep = "")
cat("  components on a DIFFERENT wild-type residue                : ",
    length(recs) - same, "\n", sep = "")
cat("written to ", out, "\n", sep = "")
