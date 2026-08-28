#### apply_lofreq_to_consensus.R
####
#### Builds a FASTA carrying a sample's low-frequency LoFreq variants, for
#### FluMut to screen. FluMut consumes sequence, not VCF, so the calls have to
#### be painted onto something.
####
#### That something is the REFERENCE, not the IRMA consensus, and the reason is
#### a coordinate bug rather than a preference.
####
#### LoFreq calls positions against the bwa-aligned reference. IRMA builds its
#### consensus de novo. The two agree only while a sample is reference-like:
#### measured on the swine WGS run, 39 of 1,123 consensus records differ in
#### length from the reference and 20 of those do not start where the reference
#### starts. Writing a reference-coordinate call onto the consensus therefore
#### mutated the wrong base in those records, silently, with nothing but a
#### bounds check between it and a wrong answer.
####
#### Painting onto the reference removes the mismatch instead of correcting it.
#### It is sound because flumut_freq_threshold defaults to 0.01, so majority
#### variants are applied too — reference + LoFreq reconstructs the sample's own
#### background. Verified: it reproduces the IRMA consensus exactly in 1,075 of
#### 1,084 comparable records (99.2%).
####
#### Plain FLUMUT still runs on the IRMA consensus, where no translation is
#### needed and de novo assembly is the more robust choice for divergent
#### samples. Only this low-frequency step moves.
####
#### NO-COVERAGE POSITIONS ARE COUNTED AND REPORTED, NOT MASKED — 2026-08-08.
####
#### The residual this addresses is real: where a sample has no reads, a
#### reference base stands in for no data, FluMut screens it, and "no marker
#### here" gets reported on the strength of nothing at all.
####
#### **Writing N there was built, measured, and rejected.** On MC-696 (166
#### uncovered positions in NA) it suppressed 6 markers correctly and
#### MANUFACTURED 3 that were not there — NA-1:I222K, NA-1:I223K, NA-1:Q136R —
#### because a heavily masked segment perturbs FluMut's own alignment and shifts
#### which residue it reads at each numbered position. Verified on one sample in
#### isolation, so it is not a cross-sample table artefact. Masking with "-"
#### was worse: 15 markers lost and 3 different spurious ones. A single N is
#### harmless and ten are harmless; 166 are not.
####
#### Inventing a marker is a worse failure than reporting one on thin evidence,
#### so the sequence is left exactly as it was. **Do not re-implement masking
#### here** — it looks like the obvious fix and it is measurably not.
####
#### What happens instead: depth is read, the exposure is counted and logged, and
#### DEPTH_PROFILE publishes the per-position depth so FluLens can flag markers
#### resting on no coverage. That is the right place for it, because FluLens is
#### the only component that maps FluMut's own numbering (HA1-5, NA-1, ...) back
#### to reference positions — see calibrateFluMut() in its handoff. A per-marker
#### coverage check needs that mapping; this script does not have it.
####
#### Depth comes from DEPTH_PROFILE (`samtools depth -a` on the reference-aligned
#### BAM), already in the reference coordinates this script paints in. IRMA's own
#### coverage tables are NOT usable: they are in IRMA's de novo consensus
#### coordinates, and reconciling those back is precisely the coordinate problem
#### reference-painting exists to avoid.
####
#### Zero coverage and below-MIN_DEPTH are counted separately. Masking at
#### depth < 100 was never on the table: MC-696 has 86% of its positions in the
#### 1-99 band and MC-533 75%, so it would blank most thin libraries outright.
####
#### Degrades: no depth directory means no counts, and everything else unchanged.
####
#### Usage:
####   Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold> \
####       <reference.fa> <depth_dir|NULL> [--min-depth=N] <irma_fasta> <lofreq_vcf> [...]
####
#### --min-depth is the run's MIN_DEPTH (default 100). Reported against, never
#### masked at. Compared to the caller-visible column when the depth files carry
#### one, since that is what MIN_DEPTH is tested against everywhere else;
#### three-column files fall back to raw depth and the summary says so.
####
#### Output headers are already >sample_SEGMENT, ready for FluMut. Do NOT pass
#### this through rename_for_flumut.R: that takes the sample name from the
#### filename, and one combined FASTA would make every record "mutated_*".

args = commandArgs(trailingOnly = TRUE)

# --min-depth=N is pulled out before anything positional is read. A flag rather
# than a fifth positional because everything past argument 4 is a variadic
# FASTA/VCF list, so a new positional would be silently swallowed as an IRMA
# FASTA by any caller not updated in step. This was previously hardcoded to 100
# and ignored the run's MIN_DEPTH.
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

# "NULL", empty, or a directory that is not there all mean "no depth source" and
# leave the sequence unmasked. Stated once here so the per-sample lookup below
# has nothing to decide.
use_depth = !is.na(depth_dir) && nzchar(depth_dir) &&
            depth_dir != "NULL" && dir.exists(depth_dir)
if (!use_depth)
  cat("No depth source: no-coverage positions will take a reference base, as before.\n",
      "  Absence of a marker is NOT evidence of absence in this output.\n",
      sep = "", file = stderr())

# MIN_DEPTH_NOTE is set from --min-depth above and reported alongside the mask,
# so the exposure below the floor stays visible without a second threshold being
# imposed on the sequence itself.

# Returns list(raw = segment -> depths, vis = segment -> depths or NULL), or NULL.
#
# Two depths, because DEPTH_PROFILE publishes two: column 3 is every aligned
# base, column 4 what the callers can actually see. Older runs are three-column
# and have no `vis`.
read_depth = function(sample_name) {
  if (!use_depth) return(NULL)
  p = file.path(depth_dir, paste0(sample_name, ".depth"))
  if (!file.exists(p) || file.info(p)$size == 0) return(NULL)
  # Column count read from the file, not assumed: colClasses is RECYCLED when
  # shorter than the row, so a hardcoded three-element vector against a
  # four-column file works by luck rather than intent.
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

# Same rule rename_for_flumut.R uses, because this writes FluMut-ready headers
# itself. It has to: that script takes the sample name from the FILENAME, and the
# low-frequency path hands it one combined mutated.fasta, so every record in every
# sample came out named "mutated_HA", "mutated_NA" and so on — 200 records
# collapsing to 8 names with sample identity gone entirely. A marker found that way
# could never be attributed to a sample.
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

  # The consensus is read ONLY to learn which segments this sample actually
  # assembled. Its sequence is deliberately not used.
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

    # Zero-coverage positions are COUNTED, never written into the sequence.
    #
    # Writing N there was implemented and then rejected on measurement: it
    # suppressed 6 markers correctly on MC-696 and MANUFACTURED 3 that were not
    # there (NA-1:I222K, NA-1:I223K, NA-1:Q136R), because a heavily masked
    # segment perturbs FluMut's own alignment and shifts which residue it reads
    # at each numbered position. Masking with "-" instead was worse: 15 lost and
    # 3 different spurious markers. A single N is harmless and ten are harmless;
    # 166 are not. Inventing a marker is a worse failure than reporting one on
    # thin evidence, so the sequence is left alone.
    #
    # The counts go to the log and the depth files are published, so the
    # exposure is visible and FluLens can flag it against its own calibrated
    # marker coordinates — which is where a per-marker coverage check belongs,
    # because that is the only place FluMut's numbering is mapped back to
    # reference positions.
    dv = if (!is.null(depth)) depth$raw[[chrom]] else NULL
    vv = if (!is.null(depth) && !is.null(depth$vis)) depth$vis[[chrom]] else NULL
    if (!is.null(dv)) {
      n = min(length(dv), length(seq_chars))
      if (n > 0) {
        # Two questions, two depths. "Was there any read here at all" is a raw
        # question and stays on column 3. "Did a caller have enough to work
        # with" is not — MIN_DEPTH is tested against the filtered count, so
        # counting thin against column 3 understates the exposure. Falls back
        # to raw for three-column files, as this always did.
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
              if (masked > 0) sprintf(" [%d position(s) with NO coverage — markers there rest on no data]", masked) else ""),
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
  cat("  Absence of a marker at those positions is NOT evidence of absence. They are\n",
      "  reported rather than masked because writing N into the sequence measurably\n",
      "  MANUFACTURES markers — see the script header. Per-position depth is published\n",
      "  alongside the results for FluLens to flag against calibrated marker positions.\n",
      sep = "", file = stderr())
}
