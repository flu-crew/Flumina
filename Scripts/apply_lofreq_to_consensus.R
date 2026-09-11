#### apply_lofreq_to_consensus.R
####
#### Generate a FASTA that carries a sample's low-frequency LoFreq variants, for
#### FluMut screening. FluMut needs sequence input, not a VCF, so the variant
#### calls are projected onto a reference sequence.
####
#### Project the calls onto the reference, not onto the IRMA consensus, to keep
#### coordinates consistent. LoFreq calls variants against a bwa-aligned
#### reference; IRMA builds a de novo consensus. Applying reference-coordinate
#### variants to an IRMA consensus can give wrong substitutions, because of indel
#### differences and different start positions between the consensus and the
#### reference.
####
#### The default flumut_freq_threshold is 0.01, so majority variants are applied
#### too. This rebuilds the sample's consensus background while it keeps the
#### reference coordinates.
####
#### The standard FluMut path still runs on the IRMA consensus; de novo assembly
#### is best for divergent samples. Only this low-frequency step uses the
#### reference projection.
####
#### Note: uncovered positions are counted and reported, but not masked. Masking
#### uncovered regions with 'N' or '-' perturbs FluMut's alignment and both
#### suppresses real markers and creates spurious ones. Instead, depth is
#### reported via DEPTH_PROFILE (`samtools depth -a` on the reference-aligned
#### BAM), so downstream tools such as FluLens can judge marker validity by read
#### depth.
####
#### IRMA's coverage tables are not used here: they are in de novo consensus
#### coordinates, which would reintroduce the coordinate mismatch.
####
#### The script degrades gracefully: if the depth directory is absent, depth
#### counts are omitted and sequence generation still runs.
####
#### Usage:
####   Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold> \
####       <reference.fa> <depth_dir|NULL> [--min-depth=N] <irma_fasta> <lofreq_vcf> [...]
####
#### --min-depth sets the minimum depth threshold (default 100). It is for
#### reporting only and does not affect sequence masking. It uses the
#### caller-visible depth column when available, or falls back to raw depth.
####
#### Output FASTA headers use the form >sample_SEGMENT for FluMut. Do not process
#### this file with rename_for_flumut.R, which derives the sample name from the
#### file name and would give ambiguous headers for a combined FASTA.

args = commandArgs(trailingOnly = TRUE)

# Extract `--min-depth=N` before reading positional arguments. It is a flag, not
# a fifth positional argument, because every argument after argument 4 forms a
# variadic FASTA/VCF list; a new positional argument would be consumed as an
# IRMA FASTA.
min_depth_flag = grep("^--min-depth=", args, value = TRUE)
MIN_DEPTH_NOTE = if (length(min_depth_flag))
  suppressWarnings(as.numeric(sub("^--min-depth=", "",
                                  min_depth_flag[length(min_depth_flag)]))) else 100
if (is.na(MIN_DEPTH_NOTE) || MIN_DEPTH_NOTE < 0) {
  cat("Error: --min-depth must be a non-negative number\n", file = stderr())
  quit(status = 1)
}
args = args[!grepl("^--min-depth=", args)]

if (length(args) < 6) {
  cat("Usage: Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold>",
      "<reference.fa> <depth_dir|NULL> [--min-depth=N] <irma_fasta> <lofreq_vcf> [...]\n",
      file = stderr())
  quit(status = 1)
}

output_fasta   = args[1]
freq_threshold = as.numeric(args[2])
reference_path = args[3]
depth_dir      = args[4]

if (is.na(freq_threshold) || freq_threshold < 0 || freq_threshold > 1) {
  cat("Error: freq_threshold must be a number between 0 and 1\n", file = stderr())
  quit(status = 1)
}
if (!file.exists(reference_path)) {
  cat("Error: reference not found:", reference_path, "\n", file = stderr())
  quit(status = 1)
}

pairs = args[-c(1, 2, 3, 4)]
if (length(pairs) %% 2 != 0) {
  cat("Error: must provide paired FASTA/VCF arguments\n", file = stderr())
  quit(status = 1)
}

# "NULL", an empty value, or a nonexistent directory means that no depth source
# is available and leaves the sequence unmasked. Normalize this once so the
# per-sample lookup has no additional cases to handle.
use_depth = !is.na(depth_dir) && nzchar(depth_dir) &&
            depth_dir != "NULL" && dir.exists(depth_dir)
if (!use_depth)
  cat("No depth source: no-coverage positions will take a reference base, as before.\n",
      "  Absence of a marker is not evidence of absence in this output.\n",
      sep = "", file = stderr())

# MIN_DEPTH_NOTE is set from `--min-depth` above and reported alongside the mask,
# so the exposure below the floor stays visible without a second threshold being
# imposed on the sequence itself.

# Return `list(raw = segment -> depths, vis = segment -> depths or NULL)`, or NULL.
#
# Store two depth values because DEPTH_PROFILE publishes both: column 3 contains
# every aligned base, and column 4 contains the bases visible to callers. Older
# runs have three columns and no `vis` value.
read_depth = function(sample_name) {
  if (!use_depth) return(NULL)
  p = file.path(depth_dir, paste0(sample_name, ".depth"))
  if (!file.exists(p) || file.info(p)$size == 0) return(NULL)
  # Read the column count instead of assuming it. `colClasses` is recycled when
  # shorter than a row, so a hard-coded three-element vector for a four-column
  # file would work accidentally rather than by design.
  first = readLines(p, n = 1, warn = FALSE)
  if (!length(first) || !nzchar(first)) return(NULL)
  ncol_d = length(strsplit(first, "\t")[[1]])
  classes = c("character", "integer", "integer")
  if (ncol_d >= 4) classes = c(classes, "integer")
  if (ncol_d > 4)  classes = c(classes, rep("NULL", ncol_d - 4))
  d = try(utils::read.table(p, sep = "\t", header = FALSE, quote = "",
                            comment.char = "", stringsAsFactors = FALSE,
                            colClasses = classes),
          silent = TRUE)
  if (inherits(d, "try-error") || nrow(d) == 0) return(NULL)
  list(raw = split(d[[3]], d[[1]]),
       vis = if (ncol_d >= 4) split(d[[4]], d[[1]]) else NULL)
}

read_fasta = function(path) {
  lines = readLines(path, warn = FALSE)
  out = list(); nm = NULL; buf = character(0)
  flush = function() {
    if (!is.null(nm)) out[[nm]] <<- paste(buf, collapse = "")
  }
  for (ln in lines) {
    if (startsWith(ln, ">")) {
      flush()
      nm = strsplit(sub("^>", "", ln), "[ \t]")[[1]][1]
      buf = character(0)
    } else if (!is.null(nm)) {
      buf = c(buf, trimws(ln))
    }
  }
  flush()
  out
}

# Use the same rule as rename_for_flumut.R, because this script writes
# FluMut-ready headers itself. That script derives the sample name from the file
# name, but the low-frequency path hands it one combined mutated.fasta, so every
# record would be named "mutated_HA", "mutated_NA", and so on, and lose its
# sample identity.
normalise_segment <- function(raw) {
  pattern = "^(?:[A-Za-z]_)?([A-Za-z0-9]+?)(?:_[HN]\\d+)?$"
  m = regexpr(pattern, raw, perl = TRUE)
  if (m == -1) return(toupper(raw))
  toupper(substr(raw, attr(m, "capture.start")[1],
                 attr(m, "capture.start")[1] + attr(m, "capture.length")[1] - 1))
}

reference = read_fasta(reference_path)
cat(sprintf("Reference: %s (%d segments)\n", reference_path, length(reference)), file = stderr())

n_samples = length(pairs) / 2
out_fh = file(output_fasta, "w")
total_applied = 0; total_skipped = 0; total_oob = 0
total_masked = 0; total_thin = 0; samples_masked = 0
# Whether any depth file carried the caller-visible column, so the summary can
# say which quantity "thin" was counted against instead of leaving it ambiguous.
any_visible = FALSE

for (i in seq(1, length(pairs), by = 2)) {
  fasta_path = pairs[i]
  vcf_path   = pairs[i + 1]

  sample_name = basename(fasta_path)
  for (suffix in c(".fasta", ".fa", ".fna"))
    if (endsWith(sample_name, suffix))
      sample_name = substr(sample_name, 1, nchar(sample_name) - nchar(suffix))

  # Read the consensus only to determine which segments this sample assembled;
  # do not use its sequence.
  assembled = names(read_fasta(fasta_path))

  variants = list()
  vcf_lines = readLines(vcf_path, warn = FALSE)
  for (line in vcf_lines) {
    if (startsWith(line, "#")) next
    f = strsplit(line, "\t")[[1]]
    if (length(f) < 8) next
    af = NA
    if (grepl("AF=", f[8])) af = as.numeric(sub("AF=", "", regmatches(f[8], regexpr("AF=[^;]+", f[8]))))
    if (is.na(af) || af < freq_threshold) next
    # Substitutions only: an indel would shift every downstream coordinate and
    # there is no alignment here to shift them against.
    if (nchar(f[4]) != 1 || nchar(f[5]) != 1) next
    variants[[length(variants) + 1]] = list(chrom = f[1], pos = as.integer(f[2]), alt = f[5])
  }

  depth = read_depth(sample_name)
  if (!is.null(depth) && !is.null(depth$vis)) any_visible = TRUE

  applied = 0; oob = 0; emitted = 0; masked = 0; thin = 0
  for (chrom in names(reference)) {
    if (!(chrom %in% assembled)) { total_skipped = total_skipped + 1; next }
    seq_chars = strsplit(reference[[chrom]], "")[[1]]
    for (v in variants) {
      if (v$chrom != chrom) next
      if (v$pos >= 1 && v$pos <= length(seq_chars)) {
        seq_chars[v$pos] = v$alt; applied = applied + 1
      } else oob = oob + 1
    }

    # Count zero-coverage positions, but do not write them into the sequence.
    #
    # Writing N (or "-") at these positions was evaluated and rejected: heavy
    # masking perturbs FluMut's own alignment and shifts which residue it reads
    # at each numbered position, which both suppresses real markers and invents
    # spurious ones. Inventing a marker is a worse failure than reporting one on
    # thin evidence, so the sequence is left alone.
    #
    # The counts are logged and the depth files are published, so the
    # exposure is visible and FluLens can flag it against its own calibrated
    # marker coordinates - which is where a per-marker coverage check belongs,
    # because that is the only place FluMut's numbering is mapped back to
    # reference positions.
    dv = if (!is.null(depth)) depth$raw[[chrom]] else NULL
    vv = if (!is.null(depth) && !is.null(depth$vis)) depth$vis[[chrom]] else NULL
    if (!is.null(dv)) {
      n = min(length(dv), length(seq_chars))
      if (n > 0) {
        # Track two quantities. Column 3 answers whether any read was present;
        # column 4 answers whether a caller had sufficient usable depth. MIN_DEPTH
        # is tested against the filtered count, so counting thin positions from
        # column 3 would understate exposure. Fall back to raw depth for
        # three-column files.
        masked = masked + sum(dv[seq_len(n)] == 0)
        tv = if (!is.null(vv) && length(vv) >= n) vv else dv
        thin = thin + sum(dv[seq_len(n)] > 0 & tv[seq_len(n)] < MIN_DEPTH_NOTE)
      }
    }

    writeLines(c(paste0(">", sample_name, "_", normalise_segment(chrom)),
                 paste(seq_chars, collapse = "")), out_fh)
    emitted = emitted + 1
  }

  total_applied = total_applied + applied; total_oob = total_oob + oob
  total_masked = total_masked + masked; total_thin = total_thin + thin
  if (masked > 0) samples_masked = samples_masked + 1
  cat(sprintf("  %s: %d variant(s) applied across %d segment(s)%s%s\n",
              sample_name, applied, emitted,
              if (oob > 0) sprintf(" [%d outside the reference, dropped]", oob) else "",
              if (masked > 0) sprintf(" [%d position(s) with NO coverage - markers there rest on no data]", masked) else ""),
      file = stderr())
}

close(out_fh)

cat(sprintf("Applied %d variants (AF >= %.1f%%) from %d sample(s) onto the reference -> %s\n",
            total_applied, freq_threshold * 100, n_samples, output_fasta), file = stderr())
if (total_skipped > 0)
  cat(sprintf("  %d sample-segment(s) skipped: not assembled by IRMA, so not emitted as bare reference\n",
              total_skipped), file = stderr())
if (total_oob > 0)
  cat(sprintf("  %d call(s) fell outside the reference and were dropped\n", total_oob), file = stderr())
if (use_depth) {
  cat(sprintf("  %d position(s) across %d sample(s) have ZERO coverage: the reference base stands in\n",
              total_masked, samples_masked), file = stderr())
  cat(sprintf("  %d further position(s) have coverage below %g, %s\n",
              total_thin, MIN_DEPTH_NOTE,
              if (any_visible)
                "measured as the callers see it (overlapping mates zeroed, quality floor applied)"
              else
                "measured as RAW depth: these depth files predate the caller-visible column, so the true exposure is larger"),
      file = stderr())
  cat("  Absence of a marker at those positions is not evidence of absence. They are\n",
      "  reported rather than masked because writing N into the sequence measurably\n",
      "  MANUFACTURES markers - see the script header. Per-position depth is published\n",
      "  alongside the results for FluLens to flag against calibrated marker positions.\n",
      sep = "", file = stderr())
}
