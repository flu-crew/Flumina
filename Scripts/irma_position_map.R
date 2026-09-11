#!/usr/bin/env Rscript

#### irma_position_map.R
####
#### Place IRMA's consensus onto this run's reference coordinates, so downstream
#### readers display it by lookup rather than re-derive it or assume a frame.
####
#### IRMA produces the consensus that FluMut screens. Painting calls above 50%
#### onto the reference is a second, independent derivation of the same quantity.
#### The two are not interchangeable: IRMA maps to its own iteratively refined
#### contig and recruits reads that BWA soft-clips at the segment termini, so it
#### carries consensus changes the reference-based callers never see.
####
#### Reading IRMA's consensus needs a coordinate frame, and a contig does not
#### automatically share the reference's. A contig can match the reference length,
#### be shorter (and shorter by a non-multiple of 3, which reads the tail of the
#### product out of frame), or be longer (an inserted base shifts every
#### downstream codon). The frame is established by alignment, per sample and per
#### segment, and computed once here.
####
#### Depth comes from IRMA's own coverage table, not the BWA pileup. min_depth is
#### applied to the alignment that PRODUCED the call. `<locus>-coverage.txt` is in
#### the contig's own coordinates and its Consensus column reconstructs the
#### contig, so the join is direct. Nothing is masked: the floor is recorded here
#### and left for the reader to apply.
####
#### Alignment note: this uses a semi-global alignment with free reference end
#### gaps (Biostrings "global-local"), so a truncated contig sits inside the
#### reference without paying for the flanks.
####
#### Outputs (see the column headers below): irma_position_map.tsv (one row per
#### sample x product frame), irma_consensus_aa.tsv (one row per codon), and
#### irma_variants.tsv (IRMA's minority calls, placed against the reference base).
####
#### usage:
####   Rscript irma_position_map.R --reference ref.fa --gtf reference_gtf \
####     --contigs IRMA-consensus-contigs --irma IRMA_results --min-depth 100 \
####     --out irma_position_map.tsv --out-aa irma_consensus_aa.tsv \
####     --out-var irma_variants.tsv

suppressMessages(library(Biostrings))

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

# Match/mismatch matrix over the DNA IUPAC alphabet, so alignment scores base
# equality just as the original aligner did (match +1, mismatch -2).
NUC.LETTERS = strsplit("ACGTRYSWKMBDHVN", "")[[1]]
NUC.MAT = matrix(-2, length(NUC.LETTERS), length(NUC.LETTERS),
                 dimnames = list(NUC.LETTERS, NUC.LETTERS))
diag(NUC.MAT) = 1

translate.codon = function(nt) {
  aa = unname(CODONS[toupper(nt)])
  if (is.na(aa)) "X" else aa
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

resolve.locus = function(name, seqs) {
  if (name %in% names(seqs)) return(name)
  for (k in names(seqs)) if (short.locus(k) == short.locus(name)) return(k)
  NULL
}

# IRMA's minority-allele calls for one segment, in CONTIG coordinates. Returns a
# data frame with one row per call, or an empty frame when the table is absent or
# missing a needed column.
read.variants = function(path) {
  empty = data.frame(position = integer(0), cons = character(0),
                     cfreq = numeric(0), minor = character(0),
                     mfreq = numeric(0), mcount = integer(0),
                     total = integer(0), stringsAsFactors = FALSE)
  if (!file.exists(path)) return(empty)
  lines = readLines(path)
  if (length(lines) < 1L) return(empty)
  head = strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  need = c("Position", "Consensus_Allele", "Minority_Allele",
           "Consensus_Frequency", "Minority_Frequency", "Minority_Count", "Total")
  if (!all(need %in% head)) return(empty)
  ix = setNames(seq_along(head), head)
  rows = list()
  for (line in lines[-1]) {
    c = strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(c) < length(head)) next
    p = suppressWarnings(as.integer(c[ix[["Position"]]]))
    cf = suppressWarnings(as.numeric(c[ix[["Consensus_Frequency"]]]))
    mf = suppressWarnings(as.numeric(c[ix[["Minority_Frequency"]]]))
    mc = suppressWarnings(as.integer(c[ix[["Minority_Count"]]]))
    tot = suppressWarnings(as.integer(c[ix[["Total"]]]))
    if (any(is.na(c(p, cf, mf, mc, tot)))) next
    rows[[length(rows) + 1L]] = data.frame(
      position = p, cons = toupper(c[ix[["Consensus_Allele"]]]), cfreq = cf,
      minor = toupper(c[ix[["Minority_Allele"]]]), mfreq = mf,
      mcount = mc, total = tot, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

# IRMA's per-position depth, by CONTIG coordinate (1-based). NULL when absent.
read.depth = function(path) {
  if (!file.exists(path)) return(NULL)
  lines = readLines(path)
  if (length(lines) < 1L) return(integer(0))
  depth = integer(0)
  for (line in lines[-1]) {
    c = strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(c) < 3L) next
    v = suppressWarnings(as.integer(c[3]))
    depth = c(depth, if (is.na(v)) 0L else v)
  }
  depth
}

colinear.identity = function(ref.seq, irma.seq) {
  same = sum(charToRaw(toupper(ref.seq)) == charToRaw(toupper(irma.seq)))
  same / max(nchar(ref.seq), 1L)
}

# Reference 0-based index -> contig 0-based index, for positions that align.
# Equal length and high identity is taken as colinear without aligning. Otherwise
# a semi-global alignment (free reference end gaps) places the contig.
map.positions = function(ref.seq, irma.seq) {
  nref = nchar(ref.seq)
  if (nref == nchar(irma.seq) && colinear.identity(ref.seq, irma.seq) >= 0.95) {
    return(list(map = seq.int(0L, nref - 1L), status = "identity"))
  }
  a = pairwiseAlignment(pattern = DNAString(toupper(irma.seq)),
                        subject = DNAString(toupper(ref.seq)),
                        type = "global-local", substitutionMatrix = NUC.MAT,
                        gapOpening = 9.5, gapExtension = 0.5)
  patc = strsplit(as.character(alignedPattern(a)), "", fixed = TRUE)[[1]]
  subc = strsplit(as.character(alignedSubject(a)), "", fixed = TRUE)[[1]]
  map = rep(NA_integer_, nref)
  ri = start(subject(a)) - 1L
  ci = start(pattern(a)) - 1L
  for (j in seq_along(subc)) {
    sg = subc[j] == "-"
    pg = patc[j] == "-"
    if (!sg && !pg) map[ri + 1L] = ci
    if (!sg) ri = ri + 1L
    if (!pg) ci = ci + 1L
  }
  list(map = map, status = "aligned")
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

reference.path = get.flag("--reference")
gtf.dir = get.flag("--gtf")
contigs.dir = get.flag("--contigs")
irma.dir = get.flag("--irma")
min.depth = suppressWarnings(as.integer(get.flag("--min-depth", "100")))
out.path = get.flag("--out", "irma_position_map.tsv")
out.aa.path = get.flag("--out-aa", "irma_consensus_aa.tsv")
out.var.path = get.flag("--out-var", "irma_variants.tsv")

if (is.null(reference.path) || is.null(gtf.dir) || is.null(contigs.dir)) {
  cat("Usage: Rscript irma_position_map.R --reference <fa> --gtf <dir> ",
      "--contigs <dir> [--irma <dir>] [--min-depth N] [--out ...]\n",
      sep = "", file = stderr())
  quit(status = 1)
}

if (!dir.exists(contigs.dir)) {
  cat("irma_position_map: no IRMA consensus contigs; nothing written\n",
      file = stderr())
  quit(status = 0)
}

seqs = read.fasta(reference.path)
cds = read.cds(gtf.dir)
if (length(cds) == 0L || length(seqs) == 0L) {
  cat("irma_position_map: no reference CDS could be read; nothing written\n",
      file = stderr())
  quit(status = 0)
}

# Flatten each product's reference CDS into 0-based reference indices, so a codon
# is three lookups rather than an interval walk.
products = list()
for (product in sort(names(cds))) {
  rec = cds[[product]]
  locus = resolve.locus(rec$seqname, seqs)
  if (is.null(locus)) {
    products[product] = list(NULL)   # keep the key with a NULL value
    next
  }
  idx = integer(0)
  for (e in rec$exons) idx = c(idx, seq.int(e[1] - 1L, e[2] - 1L))
  products[[product]] = list(locus = locus, idx = idx)
}

# Row accumulators, collected per sample and combined at the end.
map.blocks = list()
aa.blocks = list()
var.blocks = list()
tally = c(identity = 0L, aligned = 0L, absent = 0L, no_ref = 0L)
vtally = c(rows = 0L, unplaced = 0L, consensus_differs_from_reference = 0L,
           minority_is_reference = 0L)
aa.status.tally = integer(0)
bump = function(v, name, by = 1L) {
  if (is.na(v[name])) v[name] = 0L
  v[name] = v[name] + by
  v
}

fmt.g = function(x) sprintf("%.6g", x)

contig.files = sort(list.files(contigs.dir, pattern = "\\.fasta$"))
for (fn in contig.files) {
  sample = sub("\\.fasta$", "", fn)
  contigs = read.fasta(file.path(contigs.dir, fn))
  map.s = character(0)
  aa.s = character(0)
  var.s = character(0)

  # One alignment per segment, shared by every product on it (PA and PA-X share a
  # contig, as do NS1/NEP, M1/M2 and PB1/PB1-F2).
  frames = list()
  depths = list()
  for (product in names(products)) {
    rec = products[[product]]
    if (is.null(rec)) next
    locus = rec$locus
    if (locus %in% names(frames)) next
    con = get.seq(contigs, locus)
    if (is.null(con)) con = get.seq(contigs, resolve.locus(locus, contigs))
    if (is.null(con)) {
      frames[[locus]] = list(map = NULL, status = "absent")
      next
    }
    frames[[locus]] = map.positions(get.seq(seqs, locus), con)
    if (!is.null(irma.dir)) {
      depths[[locus]] = read.depth(file.path(irma.dir, sample, "tables",
                                             paste0(locus, "-coverage.txt")))
    }
  }

  # Place IRMA's minority calls on the reference and state them against the
  # reference base, using the same alignment the consensus used.
  if (!is.null(irma.dir)) {
    for (locus in sort(names(frames))) {
      fr = frames[[locus]]
      if (is.null(fr$map)) next
      # contig index -> reference index (invert the forward map).
      fwd = fr$map
      back = rep(NA_integer_, length(fwd))   # contig idx (0-based) -> ref idx
      placed = which(!is.na(fwd))
      back.len = if (length(placed)) max(fwd[placed]) + 1L else 0L
      back = rep(NA_integer_, max(back.len, 1L))
      for (r0 in placed) back[fwd[r0] + 1L] = r0 - 1L
      ref.seq = get.seq(seqs, locus)
      vars = read.variants(file.path(irma.dir, sample, "tables",
                                     paste0(locus, "-variants.txt")))
      if (nrow(vars) > 0L) for (vi in seq_len(nrow(vars))) {
        vtally = bump(vtally, "rows")
        p = vars$position[vi]
        r0 = if (p <= length(back)) back[p] else NA_integer_
        if (is.na(r0)) {
          vtally = bump(vtally, "unplaced")
          next
        }
        ref.base = toupper(substr(ref.seq, r0 + 1L, r0 + 1L))
        cons = vars$cons[vi]
        minor = vars$minor[vi]
        if (cons != ref.base) vtally = bump(vtally, "consensus_differs_from_reference")
        if (minor == ref.base) vtally = bump(vtally, "minority_is_reference")
        var.s = c(var.s, paste(sample, locus, r0 + 1L, ref.base,
          cons, fmt.g(vars$cfreq[vi]), minor, fmt.g(vars$mfreq[vi]),
          vars$mcount[vi], vars$total[vi],
          if (cons == ref.base) "yes" else "no",
          if (minor == ref.base) "yes" else "no", sep = "\t"))
      }
    }
  }

  for (product in sort(names(products))) {
    rec = products[[product]]
    if (is.null(rec)) {
      map.s = c(map.s, paste(sample, product, "", "no_ref", 0, 0, 0, 0,
                                     0, 0, min.depth, sep = "\t"))
      tally = bump(tally, "no_ref")
      next
    }

    locus = rec$locus
    fr = frames[[locus]]
    pos = fr$map
    status = fr$status
    ref.idx = rec$idx
    n.codons = length(ref.idx) %/% 3L
    ref.seq = get.seq(seqs, locus)
    con = get.seq(contigs, locus)
    if (is.null(con)) con = ""
    dep = depths[[locus]]

    if (is.null(pos)) {
      map.s = c(map.s, paste(sample, product, locus, "absent", n.codons,
                                     0, 0, 0, n.codons, 0, min.depth, sep = "\t"))
      tally = bump(tally, "absent")
      for (k in seq_len(n.codons)) {
        tri = ref.idx[(3L * (k - 1L) + 1L):(3L * (k - 1L) + 3L)]
        if (any(tri >= nchar(ref.seq))) {
          ref.aa = "?"
        } else {
          ref.aa = translate.codon(paste0(substring(ref.seq, tri + 1L, tri + 1L),
                                          collapse = ""))
        }
        aa.s = c(aa.s, paste(sample, product, k, ref.aa, "", "uncovered",
                                     "", sep = "\t"))
        aa.status.tally = bump(aa.status.tally, "uncovered")
      }
      next
    }

    counts = c(placed = 0L, changed = 0L, ambiguous = 0L, uncovered = 0L, thin = 0L)
    for (k in seq_len(n.codons)) {
      tri = ref.idx[(3L * (k - 1L) + 1L):(3L * (k - 1L) + 3L)]
      ref.aa = translate.codon(paste0(substring(ref.seq, tri + 1L, tri + 1L),
                                      collapse = ""))
      mapped = pos[tri + 1L]

      if (any(is.na(mapped)) || any(mapped >= nchar(con))) {
        counts = bump(counts, "uncovered")
        aa.s = c(aa.s, paste(sample, product, k, ref.aa, "", "uncovered",
                                     "", sep = "\t"))
        aa.status.tally = bump(aa.status.tally, "uncovered")
        next
      }

      counts = bump(counts, "placed")
      irma.aa = translate.codon(paste0(substring(con, mapped + 1L, mapped + 1L),
                                       collapse = ""))

      d = ""
      thin = FALSE
      if (!is.null(dep) && length(dep) > 0L) {
        vals = dep[mapped[mapped < length(dep)] + 1L]
        if (length(vals) > 0L) {
          d = min(vals)
          thin = d < min.depth
        }
      }

      if (irma.aa == "X") {
        counts = bump(counts, "ambiguous")
        aa.s = c(aa.s, paste(sample, product, k, ref.aa, "X",
                                     "ambiguous", d, sep = "\t"))
        aa.status.tally = bump(aa.status.tally, "ambiguous")
      } else if (irma.aa != ref.aa) {
        counts = bump(counts, "changed")
        if (thin) {
          counts = bump(counts, "thin")
          state = "change_thin"
        } else {
          state = if (!identical(d, "")) "change" else "change_nodepth"
        }
        aa.s = c(aa.s, paste(sample, product, k, ref.aa, irma.aa, state,
                                     d, sep = "\t"))
        aa.status.tally = bump(aa.status.tally, state)
      } else if (thin) {
        counts = bump(counts, "thin")
        aa.s = c(aa.s, paste(sample, product, k, ref.aa, irma.aa, "thin",
                                     d, sep = "\t"))
        aa.status.tally = bump(aa.status.tally, "thin")
      }
    }

    map.s = c(map.s, paste(sample, product, locus, status, n.codons,
      counts[["placed"]], counts[["changed"]], counts[["ambiguous"]],
      counts[["uncovered"]], counts[["thin"]], min.depth, sep = "\t"))
    tally = bump(tally, status)
  }

  map.blocks[[fn]] = map.s
  aa.blocks[[fn]] = aa.s
  var.blocks[[fn]] = var.s
}

map.lines = unlist(map.blocks, use.names = FALSE)
aa.lines = unlist(aa.blocks, use.names = FALSE)
var.lines = unlist(var.blocks, use.names = FALSE)

writeLines(c(paste("sample", "product", "locus", "status", "ref_aa_len",
                   "placed", "changed", "ambiguous", "uncovered", "thin",
                   "min_depth", sep = "\t"), map.lines), out.path)
writeLines(c(paste("sample", "product", "ref_pos", "ref_aa", "irma_aa", "status",
                   "irma_depth", sep = "\t"), aa.lines), out.aa.path)
writeLines(c(paste("sample", "locus", "ref_pos", "ref_base", "cons_allele",
                   "cons_freq", "minor_allele", "minor_freq", "minor_count",
                   "total", "cons_is_ref", "minor_is_ref", sep = "\t"),
            var.lines), out.var.path)

g = function(name) { v = aa.status.tally[name]; if (is.na(v)) 0L else v }
cat(sprintf(paste0("irma_position_map: %d samples, %d sample x product frames ",
                   "(%d colinear, %d aligned, %d absent, %d no reference)\n"),
            length(contig.files), length(map.lines), tally[["identity"]],
            tally[["aligned"]], tally[["absent"]], tally[["no_ref"]]),
    file = stderr())
if (g("change_nodepth") > 0L) {
  cat(sprintf(paste0("irma_position_map: NO IRMA coverage tables (--irma not ",
                     "given or unreadable), so the depth floor was not applied ",
                     "-- %d change_nodepth, %d ambiguous, %d uncovered\n"),
              g("change_nodepth"), g("ambiguous"), g("uncovered")), file = stderr())
} else {
  cat(sprintf(paste0("irma_position_map: min_depth %d against IRMA coverage -- ",
                     "%d change, %d change_thin, %d thin, %d ambiguous, ",
                     "%d uncovered\n"),
              min.depth, g("change"), g("change_thin"), g("thin"),
              g("ambiguous"), g("uncovered")), file = stderr())
}
if (!is.null(irma.dir)) {
  cat(sprintf(paste0("irma_position_map: %d IRMA minority calls placed ",
                     "(%d unplaced) -- %d where IRMA's consensus differs from ",
                     "the reference, %d where the MINORITY allele IS the ",
                     "reference\n"),
              length(var.lines), vtally[["unplaced"]],
              vtally[["consensus_differs_from_reference"]],
              vtally[["minority_is_reference"]]), file = stderr())
}
quit(status = 0)
