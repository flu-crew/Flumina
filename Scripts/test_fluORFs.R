#### test_fluORFs.R
####
#### Regression test for the influenza ORF definitions in fluORFs.R.
####
####   Rscript Scripts/test_fluORFs.R reference.fa [reference_gtf_dir]
####
#### With no arguments it uses the reference shipped with the repo.
####
#### Two things are checked:
####
####  1. BIOLOGY. Every product's spliced coding sequence must begin with ATG,
####     have a length divisible by three, contain no internal stop, and be
####     followed immediately by a stop codon. A wrong splice junction or a
####     wrong frame offset breaks at least one of those, which is what makes
####     this worth running — the failure mode being guarded against is a
####     coordinate that is plausible but off by one or two bases.
####
####  2. AGREEMENT WITH makeGTF.R. If a reference_gtf directory is supplied,
####     every CDS interval fluORFs.R computes must match the ones makeGTF.R
####     wrote. The two encode the same coordinates separately, so this is the
####     test that catches them drifting apart.

args <- commandArgs(trailingOnly = TRUE)

script.dir <- dirname(sub("^--file=", "",
                          grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(script.dir) || !nzchar(script.dir)) script.dir <- "."
source(file.path(script.dir, "fluORFs.R"))

ref.path <- if (length(args) >= 1) args[1] else file.path(script.dir, "..", "reference.fa")
gtf.dir  <- if (length(args) >= 2) args[2] else NA

if (!file.exists(ref.path)) stop("Reference FASTA not found: ", ref.path)
ref <- flu_read_fasta(ref.path)
cat("Reference:", ref.path, "-", length(ref), "segments\n\n")

pass <- 0L; fail <- 0L
ok <- function(label, cond) {
  good <- isTRUE(cond)
  if (good) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("  %-4s %s\n", if (good) "PASS" else "FAIL", label))
}

STOPS <- c("TAA", "TAG", "TGA")

#############################################
#### 1. Biological sanity of every CDS
#############################################
cat("--- coding sequences ---\n")
for (nm in names(ref)) {
  seq.str <- ref[[nm]]
  seq.len <- nchar(seq.str)
  for (o in flu_orfs_for(nm, seq.str)) {
    cds   <- flu_cds_seq(o, seq.str)
    label <- sprintf("%s / %-7s %4d aa", nm, o$gene, nchar(cds) %/% 3)

    ok(paste(label, "| starts with ATG"), substr(cds, 1, 3) == "ATG")
    ok(paste(label, "| length divisible by 3"), nchar(cds) %% 3 == 0)

    codons <- substring(cds, seq(1, nchar(cds) - 2, by = 3),
                             seq(3, nchar(cds),     by = 3))
    ok(paste(label, "| no internal stop codon"),
       !any(codons[-length(codons)] %in% STOPS))

    # The base immediately after the CDS should open a stop codon. Skipped when
    # the product runs to the very end of the segment.
    last.end <- o$exons$end[nrow(o$exons)]
    if (last.end + 3 <= seq.len) {
      nxt <- substr(seq.str, last.end + 1, last.end + 3)
      ok(sprintf("%s | followed by a stop codon (%s)", label, nxt), nxt %in% STOPS)
    }
  }
}

#############################################
#### 2. Position mapping invariants
#############################################
cat("\n--- position mapping ---\n")
for (nm in names(ref)) {
  seq.len <- nchar(ref[[nm]])
  for (o in flu_orfs_for(nm, ref[[nm]])) {
    lab <- sprintf("%s / %s", nm, o$gene)

    # Every CDS base must map to exactly one codon, monotonically increasing.
    all.pos <- unlist(mapply(function(s, e) s:e, o$exons$start, o$exons$end,
                             SIMPLIFY = FALSE))
    m <- flu_map_positions(o, all.pos)
    ok(paste(lab, "| every CDS base maps"), !any(is.na(m$aa_position)))
    ok(paste(lab, "| cds_position is 1..n with no gaps"),
       identical(sort(m$cds_position), seq_len(flu_cds_length(o))))
    ok(paste(lab, "| codons run 1..n/3"),
       max(m$aa_position) == flu_cds_length(o) %/% 3)

    # A single-exon product starting at nucleotide 1 is the ONLY case where
    # ceiling(position/3) is right; assert that so the distinction stays visible.
    if (nrow(o$exons) == 1 && o$exons$start[1] == 1) {
      ok(paste(lab, "| agrees with ceiling(pos/3)"),
         all(m$aa_position == ceiling(all.pos / 3)))
    } else {
      ok(paste(lab, "| does NOT agree with ceiling(pos/3)"),
         !all(m$aa_position == ceiling(all.pos / 3)))
    }

    # Positions outside the CDS must be NA, never rounded into a real codon.
    outside <- setdiff(seq_len(seq.len), all.pos)
    if (length(outside))
      ok(paste(lab, "| non-coding positions map to NA"),
         all(is.na(flu_map_positions(o, outside)$aa_position)))
  }
}

#############################################
#### 3. Agreement with makeGTF.R, when available
#############################################
if (!is.na(gtf.dir) && dir.exists(gtf.dir)) {
  cat("\n--- agreement with makeGTF.R output in", gtf.dir, "---\n")
  for (nm in names(ref)) {
    gp <- file.path(gtf.dir, paste0(nm, ".gtf"))
    if (!file.exists(gp)) { cat("  (no GTF for", nm, ")\n"); next }

    mine <- do.call(rbind, lapply(flu_orfs_for(nm, ref[[nm]]), function(o)
      data.frame(gene = o$gene, start = o$exons$start, end = o$exons$end,
                 stringsAsFactors = FALSE)))
    mine <- mine[order(mine$gene, mine$start), ]

    g <- read.table(gp, sep = "\t", comment.char = "#", quote = "",
                    stringsAsFactors = FALSE)
    g <- g[g$V3 == "CDS", ]
    theirs <- data.frame(gene  = sub('.*gene_id "([^"]+)".*', "\\1", g$V9),
                         start = g$V4, end = g$V5, stringsAsFactors = FALSE)
    theirs <- theirs[order(theirs$gene, theirs$start), ]

    ok(paste(nm, "| CDS intervals match makeGTF.R"),
       identical(paste(mine$gene, mine$start, mine$end),
                 paste(theirs$gene, theirs$start, theirs$end)))
  }
} else {
  cat("\n(no reference_gtf directory given; skipping makeGTF.R agreement check)\n")
}

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(save = "no", status = 1)
