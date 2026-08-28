#### runSNPGenie.R
#### Runs SNPGenie on per-sample LoFreq VCF files produced by the Flumina
#### variant-calling pipeline, using the GTF and FASTA files created by
#### makeGTF.R.  Per-sample results are collected and written as combined
#### summary tables in the snpGenie_results directory.

args = commandArgs(trailingOnly = TRUE)
#args = "config.cfg"

# Parse configuration file
lines <- readLines(args)
config <- list()
for (line in lines) {
  line <- trimws(line)
  if (nchar(line) != 0 && !startsWith(line, "#")) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key   <- trimws(parts[1])
      value <- trimws(parts[2])
      config[[key]] <- value
    }
  }
}

work.dir      = gsub("\"", "", config$OUTPUT_DIRECTORY)
threads       = as.numeric(gsub("\"", "", config$THREADS))
snpgenie.cfg  = gsub("\"", "", config$SNPGENIE_PATH)
if (length(snpgenie.cfg) == 0L || snpgenie.cfg == "NULL") snpgenie.cfg = NULL

# Optional minimum minor-allele frequency (decimal 0-1) passed to SNPGenie's
# --minfreq. If unset/NULL/empty, SNPGenie includes all variants (default).
min.freq = gsub("\"", "", config$MIN_ALLELE_FREQ)
if (length(min.freq) == 0L || min.freq == "NULL" || min.freq == "") {
  minfreq.arg = ""
} else {
  minfreq.arg = paste0(" --minfreq=", as.numeric(min.freq))
}

# Optional minimum read depth (coverage). SNPGenie has no --mincov option, so
# variant records with DP below this are removed from each VCF before SNPGenie
# reads it. If unset/NULL/empty, no depth filter is applied.
min.cov = gsub("\"", "", config$MIN_COVERAGE)
if (length(min.cov) == 0L || min.cov == "NULL" || min.cov == "") {
  min.coverage = 0
} else {
  min.coverage = as.numeric(min.cov)
}

vcf.dir       = paste0(work.dir, "/vcf_files")
reference.dir = paste0(work.dir, "/reference_gtf")
output.dir    = paste0(work.dir, "/snpGenie_results")
samples.dir   = paste0(output.dir, "/sample_results")

dir.create(output.dir,  showWarnings = FALSE)
dir.create(samples.dir, showWarnings = FALSE)

#############################################
#### Locate SNPGenie binary
#############################################

# 1. Check PATH (e.g. snpgenie conda env already activated)
snpgenie.path = Sys.which("snpgenie.pl")

# 2. Fall back to the path supplied in the config file
if (nchar(snpgenie.path) == 0 && !is.null(snpgenie.cfg)) {
  snpgenie.path = snpgenie.cfg
}

if (nchar(snpgenie.path) == 0 || !file.exists(snpgenie.path)) {
  stop("Could not find snpgenie.pl on PATH. ",
       "Add SNPGENIE_PATH=/path/to/snpgenie.pl to your config file, ",
       "or ensure the snpgenie environment is on PATH before running.")
}

cat("Using SNPGenie:", snpgenie.path, "\n")

#############################################
#### Load reference file lists from reference_gtf/
#############################################

gtf.files   = list.files(reference.dir, pattern = "\\.gtf$",          full.names = FALSE)
gtf.files   = gtf.files[gtf.files != "combined.gtf"]   # exclude combined
fasta.files = list.files(reference.dir, pattern = "\\.fasta$|\\.fa$", full.names = FALSE)

if (length(gtf.files) == 0) {
  stop("No GTF files found in ", reference.dir,
       ". Run makeGTF.R before runSNPGenie.R.")
}

# Reference names are the GTF filenames without extension (e.g. "A_NA_N2")
# and must match the #CHROM field in the VCF files.
reference.names = gsub("\\.gtf$", "", gtf.files)

#############################################
#### Find per-sample LoFreq VCF files
#############################################

sample.dirs = list.dirs(vcf.dir, recursive = FALSE, full.names = FALSE)
sample.dirs = sample.dirs[nchar(sample.dirs) > 0]

cat("Running SNPGenie on", length(sample.dirs), "sample(s) with",
    length(reference.names), "reference segment(s)\n")

#############################################
#### Check for combined-sequence reference files (multi-segment only)
#############################################

combined.gtf.path   = paste0(reference.dir, "/combined.gtf")
combined.fasta.path = paste0(reference.dir, "/combined.fasta")
offsets.path        = paste0(reference.dir, "/combined_offsets.tsv")

run.combined = file.exists(combined.gtf.path) &&
               file.exists(combined.fasta.path) &&
               file.exists(offsets.path)

if (run.combined) {
  offsets.df = read.table(offsets.path, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
  cat("Combined-sequence SNPGenie run enabled (",
      nrow(offsets.df), "segments )\n")
} else {
  offsets.df = NULL
}

#############################################
#### Run SNPGenie in parallel across samples
#############################################

require(foreach)
require(doParallel)

my.cluster = parallel::makeCluster(threads, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

foreach::foreach(i = seq_along(sample.dirs),
                 .packages = c("foreach")) %dopar% {

  sample.name = sample.dirs[i]
  vcf.file    = paste0(vcf.dir, "/", sample.name, "/lofreq-called-variants.vcf")

  if (!file.exists(vcf.file)) return(NULL)

  # Read VCF, skipping ## meta-header lines
  vcf.lines = readLines(vcf.file)
  skip.n    = sum(grepl("^##", vcf.lines))
  VCF = read.table(vcf.file, skip = skip.n, comment.char = "",
                   header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE)

  # Optional read-depth filter: drop variant records with DP below min.coverage
  # (parsed from the INFO field), so SNPGenie only sees adequately-covered sites.
  if (min.coverage > 0 && nrow(VCF) > 0) {
    dp = suppressWarnings(as.numeric(sub(";.*", "", sub(".*DP=", "", VCF$INFO))))
    VCF = VCF[!is.na(dp) & dp >= min.coverage, ]
  }

  sample.out = paste0(samples.dir, "/", sample.name)
  dir.create(sample.out, showWarnings = FALSE)

  # SNPGenie writes results relative to the working directory, so
  # we must setwd to the sample output directory and use a relative --outdir.
  setwd(sample.out)

  for (j in seq_along(reference.names)) {

    ref.name   = reference.names[j]
    gtf.match  = gtf.files[grep(ref.name,   gtf.files,   fixed = TRUE)]
    fasta.match= fasta.files[grep(ref.name, fasta.files, fixed = TRUE)]

    if (length(gtf.match) == 0 || length(fasta.match) == 0) next

    # Subset VCF rows to this reference segment
    sub.VCF = VCF[grep(ref.name, VCF$`#CHROM`), ]
    if (nrow(sub.VCF) == 0) next

    # Write subset VCF for SNPGenie (absolute path is fine for --snpreport)
    sub.vcf.path = paste0(sample.out, "/", ref.name, ".vcf")
    write.table(sub.VCF, file = sub.vcf.path,
                quote = FALSE, row.names = FALSE, sep = "\t")

    # Remove any pre-existing per-reference output dir so SNPGenie does not
    # abort with "Results directory already exists".
    ref.out = paste0(sample.out, "/", ref.name)
    if (dir.exists(ref.out)) unlink(ref.out, recursive = TRUE)

    # Run SNPGenie: --outdir must be a relative name; SNPGenie creates it
    # under the current working directory (sample.out).
    system(paste0(snpgenie.path,
                  " --vcfformat=2",
                  minfreq.arg,
                  " --snpreport=",  sub.vcf.path,
                  " --fastafile=",  reference.dir, "/", fasta.match,
                  " --gtffile=",    reference.dir, "/", gtf.match,
                  " --outdir=",     ref.name))

    # Move the subset VCF into the per-reference output directory
    system(paste0("mv ", sub.vcf.path,
                  " ", ref.out, "/", ref.name, ".vcf"))

  }#end j loop

  # ------------------------------------------------------------------
  # Combined-sequence SNPGenie run (multi-segment references only)
  # Build a single VCF with POS values offset to match the combined
  # FASTA, then run SNPGenie once against combined.fasta / combined.gtf.
  # ------------------------------------------------------------------
  if (run.combined && !is.null(offsets.df) && nrow(VCF) > 0) {

    combined.rows = data.frame()
    for (k in seq_len(nrow(offsets.df))) {
      seg.name   = offsets.df$segment[k]
      seg.offset = offsets.df$offset[k]
      sub.rows   = VCF[grep(seg.name, VCF$`#CHROM`, fixed = TRUE), ]
      if (nrow(sub.rows) > 0) {
        sub.rows$`#CHROM` = "combined"
        sub.rows$POS      = sub.rows$POS + seg.offset
        combined.rows     = rbind(combined.rows, sub.rows)
      }
    }

    if (nrow(combined.rows) > 0) {
      combined.vcf.path = paste0(sample.out, "/combined.vcf")
      write.table(combined.rows, file = combined.vcf.path,
                  quote = FALSE, row.names = FALSE, sep = "\t")

      combined.out = paste0(sample.out, "/combined")
      if (dir.exists(combined.out)) unlink(combined.out, recursive = TRUE)

      system(paste0(snpgenie.path,
                    " --vcfformat=2",
                    minfreq.arg,
                    " --snpreport=", combined.vcf.path,
                    " --fastafile=", combined.fasta.path,
                    " --gtffile=",   combined.gtf.path,
                    " --outdir=combined"))

      system(paste0("mv ", combined.vcf.path,
                    " ", combined.out, "/combined.vcf"))
    }
  }#end combined run

  return(NULL)

}#end foreach i

parallel::stopCluster(cl = my.cluster)

#############################################
#### Collect all SNPGenie result tables
#############################################

cat("Collecting SNPGenie results...\n")

sample.names = list.dirs(samples.dir, recursive = FALSE, full.names = FALSE)
sample.names = sample.names[nchar(sample.names) > 0]

all.codon.results   = data.frame()
all.site.results    = data.frame()
all.product.results = data.frame()
all.pop.results     = data.frame()

for (i in seq_along(sample.names)) {

  dataset.names = list.dirs(paste0(samples.dir, "/", sample.names[i]),
                             recursive = FALSE, full.names = FALSE)
  # Exclude the combined run — collected separately below
  dataset.names = dataset.names[nchar(dataset.names) > 0 & dataset.names != "combined"]

  for (j in seq_along(dataset.names)) {

    base.path = paste0(samples.dir, "/", sample.names[i], "/", dataset.names[j])

    codon.path = paste0(base.path, "/codon_results.txt")
    if (file.exists(codon.path)) {
      d = read.table(codon.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Dataset = dataset.names[j], Sample = sample.names[i], d)
      all.codon.results = rbind(all.codon.results, d)
    }

    site.path = paste0(base.path, "/site_results.txt")
    if (file.exists(site.path)) {
      d = read.table(site.path, header = TRUE)
      if (nrow(d) > 0) {
        d$file = basename(d$file)
        d = cbind(Dataset = dataset.names[j], Sample = sample.names[i], d)
        all.site.results = rbind(all.site.results, d)
      }
    }

    product.path = paste0(base.path, "/product_results.txt")
    if (file.exists(product.path)) {
      d = read.table(product.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Dataset = dataset.names[j], Sample = sample.names[i], d)
      all.product.results = rbind(all.product.results, d)
    }

    pop.path = paste0(base.path, "/population_summary.txt")
    if (file.exists(pop.path)) {
      d = read.table(pop.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Dataset = dataset.names[j], Sample = sample.names[i], d)
      all.pop.results = rbind(all.pop.results, d)
    }

  }#end j loop
}#end i loop

# Write per-segment summary tables
write.csv(all.codon.results,
          paste0(output.dir, "/codon_results_summary.csv"),
          row.names = FALSE)

write.csv(all.site.results,
          paste0(output.dir, "/site_results_summary.csv"),
          row.names = FALSE)

write.csv(all.product.results,
          paste0(output.dir, "/product_results_summary.csv"),
          row.names = FALSE)

write.csv(all.pop.results,
          paste0(output.dir, "/population_results_summary.csv"),
          row.names = FALSE)

#############################################
#### Collect combined-sequence SNPGenie results (multi-segment only)
#############################################

if (run.combined) {

  cat("Collecting combined-sequence SNPGenie results...\n")

  all.combined.codon   = data.frame()
  all.combined.site    = data.frame()
  all.combined.product = data.frame()
  all.combined.pop     = data.frame()

  for (i in seq_along(sample.names)) {

    base.path = paste0(samples.dir, "/", sample.names[i], "/combined")

    codon.path = paste0(base.path, "/codon_results.txt")
    if (file.exists(codon.path)) {
      d = read.table(codon.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Sample = sample.names[i], d)
      all.combined.codon = rbind(all.combined.codon, d)
    }

    site.path = paste0(base.path, "/site_results.txt")
    if (file.exists(site.path)) {
      d = read.table(site.path, header = TRUE)
      if (nrow(d) > 0) {
        d$file = basename(d$file)
        d = cbind(Sample = sample.names[i], d)
        all.combined.site = rbind(all.combined.site, d)
      }
    }

    product.path = paste0(base.path, "/product_results.txt")
    if (file.exists(product.path)) {
      d = read.table(product.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Sample = sample.names[i], d)
      all.combined.product = rbind(all.combined.product, d)
    }

    pop.path = paste0(base.path, "/population_summary.txt")
    if (file.exists(pop.path)) {
      d = read.table(pop.path, header = TRUE)
      d$file = basename(d$file)
      d = cbind(Sample = sample.names[i], d)
      all.combined.pop = rbind(all.combined.pop, d)
    }

  }#end i loop

  write.csv(all.combined.codon,
            paste0(output.dir, "/codon_results_combined_summary.csv"),
            row.names = FALSE)

  write.csv(all.combined.site,
            paste0(output.dir, "/site_results_combined_summary.csv"),
            row.names = FALSE)

  write.csv(all.combined.product,
            paste0(output.dir, "/product_results_combined_summary.csv"),
            row.names = FALSE)

  write.csv(all.combined.pop,
            paste0(output.dir, "/population_results_combined_summary.csv"),
            row.names = FALSE)

}#end combined collection

cat("SNPGenie analysis complete. Results in:", output.dir, "\n")
