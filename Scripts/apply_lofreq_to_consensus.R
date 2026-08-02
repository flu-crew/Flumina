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
#### RESIDUAL, and it is the reason for the presence gate below: where a sample
#### has no coverage, a reference base would stand in for no data. Eight of the
#### nine residual mismatches above are exactly that — LoFreq made no call,
#### mostly near the 3' end of NS. Whole segments are handled here (a segment the
#### sample never assembled is skipped, not emitted as bare reference). Within a
#### segment it is not, because that needs per-position depth and this process
#### stages no depth source. Do not read absence of a marker as evidence.
####
#### Usage:
####   Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold> \
####       <reference.fa> <irma_fasta> <lofreq_vcf> [<irma_fasta2> <vcf2> ...]
####
#### Output headers are already >sample_SEGMENT, ready for FluMut. Do NOT pass
#### this through rename_for_flumut.R: that takes the sample name from the
#### filename, and one combined FASTA would make every record "mutated_*".

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold>",
      "<reference.fa> <irma_fasta> <lofreq_vcf> [...]\n", file = stderr())
  quit(status = 1)
}

output_fasta   = args[1]
freq_threshold = as.numeric(args[2])
reference_path = args[3]

if (is.na(freq_threshold) || freq_threshold < 0 || freq_threshold > 1) {
  cat("Error: freq_threshold must be a number between 0 and 1\n", file = stderr())
  quit(status = 1)
}
if (!file.exists(reference_path)) {
  cat("Error: reference not found:", reference_path, "\n", file = stderr())
  quit(status = 1)
}

pairs = args[-c(1, 2, 3)]
if (length(pairs) %% 2 != 0) {
  cat("Error: must provide paired FASTA/VCF arguments\n", file = stderr())
  quit(status = 1)
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

  applied = 0; oob = 0; emitted = 0
  for (chrom in names(reference)) {
    if (!(chrom %in% assembled)) { total_skipped = total_skipped + 1; next }
    seq_chars = strsplit(reference[[chrom]], "")[[1]]
    for (v in variants) {
      if (v$chrom != chrom) next
      if (v$pos >= 1 && v$pos <= length(seq_chars)) {
        seq_chars[v$pos] = v$alt; applied = applied + 1
      } else oob = oob + 1
    }
    writeLines(c(paste0(">", sample_name, "_", normalise_segment(chrom)),
                 paste(seq_chars, collapse = "")), out_fh)
    emitted = emitted + 1
  }

  total_applied = total_applied + applied; total_oob = total_oob + oob
  cat(sprintf("  %s: %d variant(s) applied across %d segment(s)%s\n",
              sample_name, applied, emitted,
              if (oob > 0) sprintf(" [%d outside the reference, dropped]", oob) else ""),
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
