#### apply_lofreq_to_consensus.R
####
#### Generates a FASTA file incorporating low-frequency LoFreq variants for a sample,
#### formatted for FluMut screening. FluMut requires sequence input rather than VCF,
#### necessitating the projection of variant calls onto a reference sequence.
####
#### Variant calls are projected onto the reference sequence rather than the IRMA
#### consensus to ensure coordinate consistency. LoFreq identifies variants relative
#### to a bwa-aligned reference, whereas IRMA generates a de novo consensus. Applying
#### reference-coordinate variants to an IRMA consensus sequence can result in incorrect
#### allele substitutions due to insertion/deletion discrepancies and variable start
#### positions between the consensus and the reference.
####
#### Projecting variants onto the reference resolves this coordinate mismatch. Since the
#### default flumut_freq_threshold is 0.01, majority variants are also applied, effectively
#### reconstructing the sample's consensus background while preserving accurate coordinates.
####
#### The standard FluMut pipeline continues to operate on the IRMA consensus sequence,
#### as de novo assembly remains optimal for divergent samples lacking translational shifts.
#### Only this low-frequency variant screening step utilizes the reference projection method.
####
#### Note: Uncovered positions are quantified and reported, but are not masked.
####
#### This approach addresses the issue of false positive marker reporting. Previously, a
#### lack of sequence coverage resulted in the reference base being used, which FluMut
#### could incorrectly interpret as a variant marker. However, masking uncovered regions
#### with 'N' or '-' was empirically determined to cause greater alignment perturbations
#### in FluMut, leading to both marker suppression and the generation of spurious markers.
####
#### Therefore, the sequence is left unmasked. Instead, sequence depth is assessed and
#### reported via DEPTH_PROFILE (`samtools depth -a` on the reference-aligned BAM). This
#### allows downstream tools, such as FluLens, to evaluate marker validity based on read
#### depth, as FluLens contains the necessary logic to map FluMut's internal numbering
#### back to reference coordinates.
####
#### Note that IRMA's coverage tables are not utilized here, as they correspond to the
#### de novo consensus coordinates, which would reintroduce the coordinate mismatch issue.
####
#### The script degrades gracefully: if the depth directory is absent, depth counts are
#### omitted without impacting the core sequence generation.
####
#### Usage:
####   Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold> \
####       <reference.fa> <depth_dir|NULL> [--min-depth=N] <irma_fasta> <lofreq_vcf> [...]
####
#### --min-depth specifies the minimum depth threshold (default 100). This parameter is
#### utilized for reporting purposes only and does not influence sequence masking. It
#### evaluates the caller-visible depth column when available, or falls back to raw depth.
####
#### Output FASTA headers are formatted as >sample_SEGMENT to ensure compatibility with
#### FluMut. This file should not be processed by rename_for_flumut.R, as it relies on
#### file names for sample identification, which would result in ambiguous headers for a
#### combined FASTA.

args = commandArgs(trailingOnly = TRUE)

# Extract `--min-depth=N` before reading positional arguments. It is a flag rather
# than a fifth positional argument because all arguments after argument 4 form a
# variadic FASTA/VCF list. A new positional argument would otherwise be consumed
# as an IRMA FASTA by callers that had not been updated. The previous value was
# hard-coded to 100 and ignored the run's MIN_DEPTH.
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

# Use the same rule as rename_for_flumut.R because this script writes FluMut-ready headers
# itself. The other script derives the sample name from the filename, whereas the
# low-frequency path hands it one combined mutated.fasta, so every record in every
# sample came out named "mutated_HA", "mutated_NA" and so on - 200 records
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
    # Writing N at these positions was evaluated and rejected because it
    # suppressed 6 markers correctly on MC-696 but manufactured 3 that were not
    # there (NA-1:I222K, NA-1:I223K, NA-1:Q136R), because a heavily masked
    # segment perturbs FluMut's own alignment and shifts which residue it reads
    # at each numbered position. Masking with "-" instead was worse: 15 lost and
    # 3 different spurious markers. A single N is harmless and ten are harmless;
    # 166 are not. Inventing a marker is a worse failure than reporting one on
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
