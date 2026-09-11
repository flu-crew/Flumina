#!/usr/bin/env Rscript

#### flumut_position_map.R
####
#### Map FluMut's marker numbering onto this run's reference coordinates, so a
#### reader joins markers to variants by lookup rather than by an inferred offset.
####
#### A FluMut marker is written like `HA1-5:G224S`: protein HA1, residue 224, in
#### the numbering of FluMut's own reference. The marker string does not say where
#### residue 224 falls in the reference a run was mapped against. A constant shift
#### does not solve this. Small proteins carry too few markers to estimate a
#### shift, and one shift is correct only when the two proteins are colinear. HA
#### and NA differ by indels between subtypes, so no single offset holds across the
#### protein.
####
#### FluMut ships its reference sequences and their CDS annotations in
#### flumut_db.sqlite. The mapping is therefore read from the database and aligned,
#### not estimated. Protein-to-protein alignment handles the indels, which is why
#### HA and NA place across subtypes and not only on H5N1.
####
#### Alignment note: this uses a semi-global protein alignment (BLOSUM62) with free
#### end gaps on THIS run's ORF (Biostrings "global-local"), so a FluMut subunit
#### such as HA1 or HA2 sits inside the ORF without paying for the flanks. Internal
#### gaps are scored normally, because those are the real indels between subtypes.
####
#### Output: flumut_position_map.tsv, one row per FluMut residue that aligns.
####   label       marker prefix as it appears in results (HA1-5, NS-1, PB2)
####   protein     FluMut protein the label belongs to, per the database
####   flumut_pos  residue number as written in the marker
####   flumut_aa   residue FluMut's reference carries there
####   product     product name in THIS run's GTF
####   ref_pos     codon number in THIS run's CDS  (the join key)
####   ref_aa      residue THIS run's reference carries there
####   identical   whether the two references agree at that residue
####
#### A position absent from the file did not align and must not be placed.
####
#### usage:
####   Rscript flumut_position_map.R --reference ref.fa --gtf reference_gtf \
####     [--db flumut_db.sqlite] [--out flumut_position_map.tsv]

suppressMessages(library(Biostrings))
suppressMessages(library(DBI))

# Standard genetic code, built in the same T,C,A,G order the pipeline uses.
CODON.BASES = c("T", "C", "A", "G")
CODON.AA = strsplit(paste0("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS",
                           "RRVVVVAAAADDEEGGGG"), "")[[1]]
CODON.NAMES = character(64)
.k = 1L
for (.i in CODON.BASES) for (.j in CODON.BASES) for (.l in CODON.BASES) {
  CODON.NAMES[.k] = paste0(.i, .j, .l); .k = .k + 1L
}
CODONS = setNames(CODON.AA, CODON.NAMES)

# BLOSUM62, the same substitution matrix FluMut's own aligner uses.
data("BLOSUM62", package = "Biostrings")

# Map each FluMut protein to its product name in the GTF. FluMut splits HA into
# two mature subunits and numbers each from 1; the GTF has one HA ORF, so both
# subunits map to it and alignment determines their positions.
PRODUCT = c(HA1 = "HA", HA2 = "HA", "NA" = "NA", "NS-1" = "NS1", "NS-2" = "NEP",
            M1 = "M1", M2 = "M2", NP = "NP", PA = "PA", "PA-X" = "PA-X",
            PB1 = "PB1", "PB1-F2" = "PB1-F2", PB2 = "PB2")

# Translate a spliced CDS to protein, codon by codon, dropping a trailing partial
# codon. An unknown codon becomes X, matching FluMut's own translation.
translate.cds = function(cds) {
  cds = toupper(cds)
  n = nchar(cds)
  if (n < 3L) return("")
  starts = seq.int(1L, n - 2L, by = 3L)
  aa = vapply(starts, function(i) {
    a = unname(CODONS[substr(cds, i, i + 2L)])
    if (is.na(a)) "X" else a
  }, character(1))
  paste0(aa, collapse = "")
}

# Strip a trailing subtype suffix only: HA_H3 -> HA, NA_N2 -> NA. Leaves PA-X,
# PB1-F2, NS1, M1 and M2 untouched, matching makeGTF.R and the viewer.
norm.gene = function(g) sub("_[HN]?[0-9]+$", "", as.character(g))

short.locus = function(x) {
  x = if (length(x) == 0L || is.na(x)) "" else as.character(x)
  sub("_[A-Z][0-9]+$", "", sub("^A_", "", x))
}

# Read a FASTA into a named character vector: header first token -> sequence.
read.fasta = function(path) {
  lines = readLines(path)
  seqs = list()
  nm = NULL
  buf = character(0)
  for (line in lines) {
    line = trimws(line)
    if (!nzchar(line)) next
    if (substr(line, 1L, 1L) == ">") {
      if (!is.null(nm)) seqs[[nm]] = paste0(buf, collapse = "")
      nm = strsplit(substr(line, 2L, nchar(line)), "[[:space:]]")[[1]][1]
      buf = character(0)
    } else {
      buf = c(buf, line)
    }
  }
  if (!is.null(nm)) seqs[[nm]] = paste0(buf, collapse = "")
  unlist(seqs)
}

get.seq = function(seqs, name) {
  if (!is.null(name) && length(name) == 1L && name %in% names(seqs))
    unname(seqs[[name]]) else NULL
}

resolve.locus = function(name, seqs) {
  if (name %in% names(seqs)) return(name)
  for (k in names(seqs)) if (short.locus(k) == short.locus(name)) return(k)
  NULL
}

# CDS intervals per product, exactly as makeGTF.R defines them. combined.gtf is
# skipped and identical intervals de-duplicated, because a repeated CDS record
# would otherwise double a product's length.
read.cds = function(gtf.dir) {
  cds = list()
  seen = character(0)
  files = sort(list.files(gtf.dir, pattern = "\\.gtf$"))
  for (fn in files) {
    if (grepl("combined\\.gtf$", fn)) next
    for (line in readLines(file.path(gtf.dir, fn))) {
      if (!nzchar(trimws(line)) || substr(line, 1L, 1L) == "#") next
      c = strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(c) < 9L || c[3] != "CDS") next
      m = regmatches(c[9], regexec("gene_id \"([^\"]+)\"", c[9]))[[1]]
      if (length(m) == 0L) next
      g = norm.gene(m[2])
      key = paste(g, c[1], c[4], c[5], sep = "\t")
      if (key %in% seen) next
      seen = c(seen, key)
      if (is.null(cds[[g]])) cds[[g]] = list(seqname = c[1], exons = list())
      cds[[g]]$exons[[length(cds[[g]]$exons) + 1L]] =
        c(as.integer(c[4]), as.integer(c[5]))
    }
  }
  for (g in names(cds)) {
    ex = cds[[g]]$exons
    starts = vapply(ex, function(e) e[1], integer(1))
    cds[[g]]$exons = ex[order(starts)]
  }
  cds
}

# Translate this run's CDS exactly as makeGTF.R defines them, returning a named
# vector product -> protein.
read.our.proteins = function(reference.fa, gtf.dir) {
  cds = read.cds(gtf.dir)
  seqs = read.fasta(reference.fa)
  prot = character(0)
  for (g in names(cds)) {
    rec = cds[[g]]
    locus = resolve.locus(rec$seqname, seqs)
    if (is.null(locus)) next
    seq = get.seq(seqs, locus)
    nt = character(0)
    for (e in rec$exons) nt = c(nt, substr(seq, e[1], e[2]))
    prot[g] = translate.cds(paste0(nt, collapse = ""))
  }
  prot
}

# Splice and translate each protein from FluMut's own reference sequences.
# Returns the proteins, the marker labels (prefix -> protein), and the db version.
read.flumut.proteins = function(db) {
  con = dbConnect(RSQLite::SQLite(), dbname = db)
  on.exit(dbDisconnect(con), add = TRUE)

  ref.tab = dbGetQuery(con, 'select name, sequence from "references"')
  refs = setNames(ref.tab$sequence, ref.tab$name)

  ann.tab = dbGetQuery(con,
    "select start, end, protein_name, reference_name from annotations")
  ann = list()
  for (i in seq_len(nrow(ann.tab))) {
    p = ann.tab$protein_name[i]
    if (is.null(ann[[p]])) ann[[p]] = list(ref = ann.tab$reference_name[i],
                                           exons = list())
    ann[[p]]$exons[[length(ann[[p]]$exons) + 1L]] =
      c(ann.tab$start[i], ann.tab$end[i])
  }

  mut.tab = dbGetQuery(con, "select name, protein_name from mutations")
  labels = list()
  for (i in seq_len(nrow(mut.tab))) {
    key = strsplit(mut.tab$name[i], ":", fixed = TRUE)[[1]][1]
    if (is.null(labels[[key]])) labels[[key]] = mut.tab$protein_name[i]
  }

  version = tryCatch({
    v = dbGetQuery(con, "select major, minor from db_version")
    if (nrow(v) >= 1L) paste(v$major[1], v$minor[1], sep = ".") else NULL
  }, error = function(e) NULL)

  prot = character(0)
  for (p in names(ann)) {
    seq = if (ann[[p]]$ref %in% names(refs)) unname(refs[[ann[[p]]$ref]]) else NULL
    if (is.null(seq)) next
    ex = ann[[p]]$exons
    starts = vapply(ex, function(e) e[1], integer(1))
    ex = ex[order(starts)]
    nt = character(0)
    for (e in ex) nt = c(nt, substr(seq, e[1], e[2]))
    prot[p] = translate.cds(paste0(nt, collapse = ""))
  }
  list(prot = prot, labels = labels, version = version)
}

# Align FluMut's protein onto ours, returning a data frame of aligned index pairs
# (our_idx0, fm_idx0). FluMut's protein is the pattern and must align in full; our
# ORF is the subject and carries free end gaps ("global-local"), so a subunit sits
# inside the ORF without paying for the flanks. Gap scores reproduce the original
# aligner (open -11, extend -1): Biostrings charges gapOpening + gapExtension per
# gap position, so gapOpening 10 and gapExtension 1 give the same affine cost.
align = function(our.seq, fm.seq) {
  a = pairwiseAlignment(pattern = AAString(fm.seq), subject = AAString(our.seq),
                        type = "global-local", substitutionMatrix = BLOSUM62,
                        gapOpening = 10, gapExtension = 1)
  patc = strsplit(as.character(alignedPattern(a)), "", fixed = TRUE)[[1]]
  subc = strsplit(as.character(alignedSubject(a)), "", fixed = TRUE)[[1]]
  our.idx = integer(0)
  fm.idx = integer(0)
  ri = start(subject(a)) - 1L   # our index, 0-based
  qi = start(pattern(a)) - 1L   # fm index, 0-based
  for (j in seq_along(subc)) {
    sg = subc[j] == "-"
    pg = patc[j] == "-"
    if (!sg && !pg) {
      our.idx = c(our.idx, ri)
      fm.idx = c(fm.idx, qi)
    }
    if (!sg) ri = ri + 1L
    if (!pg) qi = qi + 1L
  }
  data.frame(our = our.idx, fm = fm.idx)
}

# Read a `--flag value` pair from the argument vector.
args = commandArgs(trailingOnly = TRUE)
get.flag = function(name, default = NULL) {
  hit = which(args == name)
  if (length(hit) == 0L) return(default)
  i = hit[length(hit)]
  if (i == length(args)) {
    cat("Error: ", name, " needs a value\n", sep = "", file = stderr())
    quit(status = 1)
  }
  args[i + 1L]
}

# Locate flumut_db.sqlite the same way FluMut does: an explicit path, else the
# copy shipped inside the flumutdb Python package on this environment's PATH.
find.db = function(explicit) {
  if (!is.null(explicit)) return(explicit)
  py = Sys.which("python3")
  if (nzchar(py)) {
    code = paste0("import flumutdb, os; ",
                  "print(os.path.join(os.path.dirname(flumutdb.__file__), ",
                  "'flumut_db.sqlite'))")
    out = suppressWarnings(system2(py, c("-c", shQuote(code)),
                                   stdout = TRUE, stderr = FALSE))
    if (length(out) >= 1L && nzchar(out[1]) && file.exists(out[1])) return(out[1])
  }
  NULL
}

reference.path = get.flag("--reference")
gtf.dir = get.flag("--gtf")
db.path = get.flag("--db")
out.path = get.flag("--out", "flumut_position_map.tsv")

if (is.null(reference.path) || is.null(gtf.dir)) {
  cat("Usage: Rscript flumut_position_map.R --reference <fa> --gtf <dir> ",
      "[--db <flumut_db.sqlite>] [--out ...]\n", sep = "", file = stderr())
  quit(status = 1)
}

db = find.db(db.path)
if (is.null(db)) {
  cat("flumut_position_map: flumut database not found; no position map written\n",
      file = stderr())
  quit(status = 0)
}

our = read.our.proteins(reference.path, gtf.dir)
fm.data = read.flumut.proteins(db)
fm = fm.data$prot
labels = fm.data$labels
dbver = fm.data$version

if (length(our) == 0L) {
  cat("flumut_position_map: no reference proteins could be translated; ",
      "no position map written\n", sep = "", file = stderr())
  quit(status = 0)
}

rows = character(0)
summary = character(0)
for (label in sort(names(labels))) {
  protein = labels[[label]]
  product = if (protein %in% names(PRODUCT)) unname(PRODUCT[[protein]]) else NULL
  fm.seq = if (protein %in% names(fm)) unname(fm[[protein]]) else NULL
  our.seq = if (!is.null(product) && product %in% names(our)) unname(our[[product]]) else NULL

  if (is.null(fm.seq) || !nzchar(fm.seq) || is.null(our.seq) || !nzchar(our.seq)) {
    msg = if (is.null(fm.seq) || !nzchar(fm.seq)) "no FluMut protein"
          else paste0("no ", product, " in this reference")
    summary = c(summary, sprintf("  %-8s %-8s -> %s", label, protein, msg))
    next
  }

  pairs = align(our.seq, fm.seq)
  our.aa = substring(our.seq, pairs$our + 1L, pairs$our + 1L)
  fm.aa = substring(fm.seq, pairs$fm + 1L, pairs$fm + 1L)
  ident = sum(our.aa == fm.aa)
  rows = c(rows, paste(label, protein, pairs$fm + 1L, fm.aa, product,
                       pairs$our + 1L, our.aa,
                       ifelse(our.aa == fm.aa, "yes", "no"), sep = "\t"))

  offsets = unique(pairs$our - pairs$fm)
  shape = if (length(offsets) == 1L) sprintf("offset %+d", offsets)
          else sprintf("%d blocks (indels)", length(offsets))
  summary = c(summary, sprintf(
    "  %-8s %-8s -> %-7s %4d/%d aligned, %5.1f%% identical, %s",
    label, protein, product, nrow(pairs), nchar(fm.seq),
    100 * ident / max(nrow(pairs), 1L), shape))
}

writeLines(c(paste("label", "protein", "flumut_pos", "flumut_aa", "product",
                   "ref_pos", "ref_aa", "identical", sep = "\t"), rows), out.path)

n.labels = if (length(rows) == 0L) 0L else length(unique(sub("\t.*$", "", rows)))
cat(sprintf("flumut_position_map: flumut db %s, %d positions across %d labels\n",
            if (is.null(dbver)) "unknown" else dbver, length(rows), n.labels),
    file = stderr())
if (length(summary) > 0L) cat(paste(summary, collapse = "\n"), "\n", sep = "",
                              file = stderr())
quit(status = 0)
