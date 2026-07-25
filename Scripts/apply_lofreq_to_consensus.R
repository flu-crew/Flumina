#### Apply LoFreq low-frequency variants to IRMA consensus sequences
#### Creates mutated FASTA for FluMut screening at user-defined variant frequency thresholds

#### Reads:
####   - IRMA consensus FASTA (one per sample)
####   - LoFreq VCF (one per sample, with AF field for allele frequency)
#### Writes:
####   - Mutated consensus FASTA containing variants above frequency threshold
####
#### Each FASTA sequence is read entirely into memory, variants are applied in-place,
#### and the mutated sequence is written to the output FASTA. Variants are applied
#### in coordinate order to ensure deterministic behavior.

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript apply_lofreq_to_consensus.R <output.fasta> <freq_threshold> ",
      "<irma_fasta> <lofreq_vcf> [<irma_fasta2> <lofreq_vcf2> ...]\n", file = stderr())
  quit(status = 1)
}

output_fasta = args[1]
freq_threshold = as.numeric(args[2])

if (is.na(freq_threshold) || freq_threshold < 0 || freq_threshold > 1) {
  cat("Error: freq_threshold must be a number between 0 and 1\n", file = stderr())
  quit(status = 1)
}

pairs = args[-c(1, 2)]
if (length(pairs) %% 2 != 0) {
  cat("Error: must provide paired FASTA/VCF arguments\n", file = stderr())
  quit(status = 1)
}

n_samples = length(pairs) / 2
out_fh = file(output_fasta, "w")
total_variants_applied = 0

for (i in seq(1, length(pairs), by = 2)) {
  fasta_path = pairs[i]
  vcf_path = pairs[i + 1]

  sample_name = basename(fasta_path)
  for (suffix in c(".fasta", ".fa", ".fna")) {
    if (endsWith(sample_name, suffix)) {
      sample_name = substr(sample_name, 1, nchar(sample_name) - nchar(suffix))
      break
    }
  }

  cat(sprintf("Processing %s: %s + %s\n", sample_name, basename(fasta_path), basename(vcf_path)),
      file = stderr())

  #############################################
  ### Parse VCF: extract variants above threshold
  #############################################

  vcf_fh = file(vcf_path, "r")
  variants = list()

  while (length(line <- readLines(vcf_fh, n = 1)) > 0) {
    if (startsWith(line, "#")) next

    fields = strsplit(line, "\t")[[1]]
    if (length(fields) < 8) next

    chrom = fields[1]
    pos = as.numeric(fields[2])
    ref = fields[4]
    alt = fields[5]
    info = fields[8]

    af = NA
    if (grepl("AF=", info)) {
      af_str = regmatches(info, regexpr("AF=[^;]+", info))
      af = as.numeric(sub("AF=", "", af_str))
    }

    if (!is.na(af) && af >= freq_threshold) {
      key = paste0(chrom, ":", pos)
      variants[[key]] = list(pos = pos, chrom = chrom, ref = ref, alt = alt, af = af)
    }
  }
  close(vcf_fh)

  n_variants_sample = length(variants)
  total_variants_applied = total_variants_applied + n_variants_sample

  if (n_variants_sample == 0) {
    cat(sprintf("  No variants above AF >= %.1f%% found\n", freq_threshold * 100),
        file = stderr())
    next
  }

  #############################################
  ### Read consensus FASTA and apply variants
  #############################################

  fasta_fh = file(fasta_path, "r")
  current_header = NULL
  current_seq = ""
  current_chrom = NULL

  while (length(line <- readLines(fasta_fh, n = 1)) > 0) {
    if (startsWith(line, ">")) {
      if (!is.null(current_seq) && nchar(current_seq) > 0) {
        seq_chars = strsplit(current_seq, "")[[1]]
        if (!is.null(current_chrom)) {
          variant_keys = names(variants)[grep(paste0("^", current_chrom, ":"), names(variants))]
          for (key in sort(variant_keys)) {
            v = variants[[key]]
            pos = v$pos
            if (pos >= 1 && pos <= length(seq_chars)) {
              seq_chars[pos] = substr(v$alt, 1, 1)
            }
          }
        }
        mutated_seq = paste(seq_chars, collapse = "")
        writeLines(c(current_header, mutated_seq), out_fh)
      }

      current_header = line
      current_chrom = strsplit(substr(line, 2, nchar(line)), " ")[[1]][1]
      current_seq = ""
    } else {
      current_seq = paste0(current_seq, line)
    }
  }

  if (!is.null(current_seq) && nchar(current_seq) > 0) {
    seq_chars = strsplit(current_seq, "")[[1]]
    if (!is.null(current_chrom)) {
      variant_keys = names(variants)[grep(paste0("^", current_chrom, ":"), names(variants))]
      for (key in sort(variant_keys)) {
        v = variants[[key]]
        pos = v$pos
        if (pos >= 1 && pos <= length(seq_chars)) {
          seq_chars[pos] = substr(v$alt, 1, 1)
        }
      }
    }
    mutated_seq = paste(seq_chars, collapse = "")
    writeLines(c(current_header, mutated_seq), out_fh)
  }

  close(fasta_fh)

  cat(sprintf("  Applied %d variant(s) to %d sequence(s)\n", n_variants_sample,
              length(grep("^>", readLines(fasta_path)))), file = stderr())
}

close(out_fh)

cat(sprintf("Applied %d total variants (AF >= %.1f%%) across %d sample(s) to %s\n",
            total_variants_applied, freq_threshold * 100, n_samples, output_fasta),
    file = stderr())
