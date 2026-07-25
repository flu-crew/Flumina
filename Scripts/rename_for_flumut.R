#### Rewrite IRMA consensus FASTA headers into sample_segment form FluMut expects
#### Concatenate all samples into one batch FASTA for FluMut screening

#### IRMA's per-sample consensus (one file per sample, e.g., <sample>.fasta)
#### carries bare headers with no sample name, e.g.:
####   >A_HA_H5
####   >A_NA_N1
####   >A_MP
####
#### FluMut expects headers like >sample_HA (sample_segment format)
#### Sample name is injected from filename; segment names are normalized
#### to canonical 8-segment codes: PB2, PB1, PA, HA, NP, NA, MP, NS

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript rename_for_flumut.R <output.fasta> <sample1.fasta> [sample2.fasta ...]\n", file = stderr())
  quit(status = 1)
}

out_path = args[1]
in_paths = args[-1]

normalise_segment <- function(raw) {
  pattern = "^(?:[A-Za-z]_)?([A-Za-z0-9]+?)(?:_[HN]\\d+)?$"
  m = regexpr(pattern, raw, perl = TRUE)
  if (m == -1) {
    return(toupper(raw))
  }
  segment = toupper(substr(raw, attr(m, "capture.start")[1],
                           attr(m, "capture.start")[1] + attr(m, "capture.length")[1] - 1))
  return(segment)
}

n_records = 0
out_fh = file(out_path, "w")

for (path in in_paths) {
  sample = basename(path)
  for (suffix in c(".fasta", ".fa", ".fna")) {
    if (endsWith(sample, suffix)) {
      sample = substr(sample, 1, nchar(sample) - nchar(suffix))
      break
    }
  }

  in_fh = file(path, "r")
  while (length(line <- readLines(in_fh, n = 1)) > 0) {
    if (startsWith(line, ">")) {
      raw_segment = strsplit(substr(line, 2, nchar(line)), " ")[[1]][1]
      segment = normalise_segment(raw_segment)
      writeLines(paste0(">", sample, "_", segment), out_fh)
      n_records = n_records + 1
    } else {
      writeLines(line, out_fh)
    }
  }
  close(in_fh)
}

close(out_fh)
cat(sprintf("wrote %d records from %d sample(s) to %s\n", n_records, length(in_paths), out_path),
    file = stderr())
