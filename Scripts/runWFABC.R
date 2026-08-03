#### runWFABC.R
#### Estimates per-site selection coefficients (s) and effective population
#### size (Ne) from allele-frequency time series using WFABC (Foll et al. 2015).
####
#### Prior to this script, run the Flumina pipeline to obtain per-sample LoFreq
#### VCF files (vcf_files/) and a metadata file that records, for every sample,
#### the individual it belongs to and the time point it was collected at.
####
#### Required R packages: data.table, Biostrings, seqinr, coda
####                      (ggplot2 optional, for the summary figure)
#### Required external tool: WFABC v1.2 (wfabc_1 and wfabc_2 binaries)
####
#### Config keys used (see config.cfg):
####   OUTPUT_DIRECTORY      - Flumina analysis directory (holds vcf_files/)
####   REFERENCE_FILE        - reference FASTA used for variant calling
####   METADATA              - CSV with a "Sample" column to join on
####   INDIVIDUAL_COLUMN     - metadata column identifying the individual that
####                           is followed through time (e.g. Animal_ID)
####   TIME_COLUMN           - metadata column with the time point (numeric)
####   GROUP_NAMES           - metadata column used to group/colour the summary
####   GENERATIONS_PER_TIME  - generations elapsed per unit of TIME_COLUMN
####   MIN_DEPTH             - minimum read depth to keep a variant
####   MIN_ALLELE_FREQUENCY  - minimum allele frequency to keep a variant
####   WFABC_PATH            - directory containing wfabc_1/wfabc_2 (optional if
####                           the binaries are already on PATH)

args = commandArgs(trailingOnly = TRUE)
#args = "config.cfg"

#############################################
#### Parse configuration file
#############################################

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

# Helper to pull a clean value from the config list
cfg = function(key, default = NULL) {
  val = config[[key]]
  if (is.null(val)) return(default)
  val = gsub("\"", "", val)
  if (length(val) == 0L || val == "" || val == "NULL") return(default)
  val
}

work.dir         = cfg("OUTPUT_DIRECTORY")
vcf.directory    = paste0(work.dir, "/vcf_files")
output.directory = paste0(work.dir, "/wfabc_analysis")

reference.path   = cfg("REFERENCE_FILE")
metadata.file    = cfg("METADATA")

# Shared influenza ORF definitions — the same ones convertVCFtoTable.R and
# findAAChanges.R use. This script carries its own copy of their VCF-to-table
# and codon-translation logic, so it needs the same correction.
wfabc.script.dir = dirname(sub("^--file=", "",
                               grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(wfabc.script.dir) || !nzchar(wfabc.script.dir)) wfabc.script.dir = "."
source(file.path(wfabc.script.dir, "fluORFs.R"))

# Metadata column names (sanitised the same way read.csv does on import)
individual.column   = make.names(cfg("INDIVIDUAL_COLUMN", "Individual"))
time.column         = make.names(cfg("TIME_COLUMN", "Time"))
group.column        = make.names(cfg("GROUP_NAMES", "Group"))
generations.per.time = as.numeric(cfg("GENERATIONS_PER_TIME", "1"))

min.depth             = as.numeric(cfg("MIN_DEPTH", "100"))
min.allele.frequency  = as.numeric(cfg("MIN_ALLELE_FREQUENCY", "0.015"))

# WFABC hangs on sites where the allele is pinned at a frequency boundary (near
# fixation/loss): the ABC sampler's acceptance rate approaches zero so wfabc_2
# can run effectively forever. Two guards handle this (see the run loop below):
#  1. Pre-filter: skip a site if the allele is near-fixed/near-lost at EVERY
#     time point (no identifiable selection; boundary-crossing swings are kept).
#  2. Timeout watchdog: kill+skip any individual wfabc call exceeding the limit.
fixation.cutoff = as.numeric(cfg("FIXATION_CUTOFF", "0.98"))  # all freq >= cutoff or <= 1-cutoff
wfabc.timeout   = as.numeric(cfg("WFABC_TIMEOUT", "120"))     # seconds per wfabc call
max.refine.iter = as.numeric(cfg("MAX_REFINE_ITER", "10"))    # cap on the ks refinement loop

# How to handle multiple variant records at the same time point ("merge" sums
# the counts; anything else keeps the records as-is and drops singletons).
duplicate.handling = "merge"

#############################################
#### Locate the WFABC binaries
#############################################

# 1. Check PATH, 2. fall back to the directory given by WFABC_PATH in config.
wfabc1.bin = Sys.which("wfabc_1")
wfabc2.bin = Sys.which("wfabc_2")
wfabc.cfg  = cfg("WFABC_PATH")
if ((nchar(wfabc1.bin) == 0 || nchar(wfabc2.bin) == 0) && !is.null(wfabc.cfg)) {
  wfabc1.bin = file.path(wfabc.cfg, "wfabc_1")
  wfabc2.bin = file.path(wfabc.cfg, "wfabc_2")
}
if (nchar(wfabc1.bin) == 0 || nchar(wfabc2.bin) == 0 ||
    !file.exists(wfabc1.bin) || !file.exists(wfabc2.bin)) {
  stop("Could not find wfabc_1/wfabc_2. Add WFABC_PATH=/path/to/wfabc/binaries ",
       "to your config file, or place the binaries on PATH before running.")
}
cat("Using WFABC:", wfabc1.bin, "/", wfabc2.bin, "\n")

dir.create(output.directory, showWarnings = FALSE)

#############################################
#### Convert LoFreq VCF files to a table
#############################################

#the string or name of the VCF file for data analysis
vcf.string = "lofreq-called-variants.vcf"

#Get multifile databases together
all.files = list.files(vcf.directory, recursive = T)
vcf.files = all.files[grep(paste0(vcf.string, "$"), all.files)]

#Collects the super cool data
header.data = c("method", "sample", "locus", "position", "reference",
                "alternative", "quality", "depth", "map_quality", "allele_frequency", "aa_position")

#Sets up data collection data.frame
collect.data = data.table::data.table(matrix(as.numeric(0),
                                             nrow = length(vcf.files)*1000,
                                             ncol = length(header.data)))
data.table::setnames(collect.data, header.data)

collect.data[, method:=as.character(method)]
collect.data[, sample:=as.character(sample)]
collect.data[, locus:=as.character(locus)]
collect.data[, reference:=as.character(reference)]
collect.data[, alternative:=as.character(alternative)]

#Loops through each locus and does operations on them
x = 1
for (i in 1:length(vcf.files)){

  #Counts comment lines to find first line
  VCF = file(paste0(vcf.directory, "/", vcf.files[i]), "r")
  skip = 0
  line = readLines(VCF, 1)

  #Finds contig line
  while(!grepl("#CHROM", line)) {
    skip = skip + 1
    line = readLines(VCF, 1)
  }

  close(VCF)

  #Reads in VCF after finding which lines to skip
  VCF = read.table(paste0(vcf.directory, "/", vcf.files[i]), skip = skip, comment.char = "", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)

  if (is.null(nrow(VCF)) == TRUE || nrow(VCF) == 0){

    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "LoFreq")
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("aa_position", header.data), value = 0)
    x = x + 1
    next
  }

  for (j in 1:nrow(VCF)){

    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "LoFreq")
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = VCF$'#CHROM'[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = VCF$POS[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = VCF$REF[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = VCF$ALT[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = VCF$QUAL[j] )

    #find depth
    depth.val = as.numeric(gsub(";", "", gsub(";.*", "", gsub(".*DP=", "", VCF[j,]$INFO))) )
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = depth.val)

    #find map quality
    if (length(grep("MQ=", VCF[j,]$INFO)) != 0){
      mq.val = as.numeric(gsub(";.*", "", gsub(".*;MQ=", "", VCF[j,]$INFO)))
    } else { mq.val = NA }
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = mq.val)

    #find allele frequency
    freq.val = as.numeric(gsub(";.*", "", gsub(".*;AF=", "", VCF[j,]$INFO)))
    data.table::set(collect.data, i = as.integer(x), j = match("allele_frequency", header.data), value = freq.val)

    # aa_position is filled in after the loop, from the real coding intervals
    # (Scripts/fluORFs.R). Unlike the main variant table this one stays at ONE
    # ROW PER CALL — the trajectories below are keyed on locus + nucleotide
    # position, and expanding a position into two products would duplicate the
    # series. Only the PRIMARY product is annotated here.
    x = x + 1
  }#end j loop

}#end i loop

#Removes empty samples
collect.data = collect.data[collect.data$sample != 0,]
collect.data = collect.data[collect.data$locus != 0,]

#Sometimes the amino acid T will turn into TRUE
collect.data$alternative[collect.data$alternative == "TRUE"] = "T"
collect.data$reference[collect.data$reference == "TRUE"] = "T"

#Annotate against the primary product's real coding interval. A call outside it
#(past the M1 or NS1 stop codon, say) gets NA rather than a fabricated codon
#number; the trajectory itself is unaffected, since that is keyed on the
#nucleotide position.
collect.data = as.data.frame(collect.data, stringsAsFactors = FALSE)
wfabc.reference = flu_read_fasta(reference.path)
collect.data = flu_annotate_positions(collect.data, wfabc.reference,
                                      locus.col = "locus", pos.col = "position",
                                      primary.only = TRUE)
cat(sprintf("WFABC amino-acid annotation: %d of %d calls fall inside their primary ORF\n",
            sum(!is.na(collect.data$aa_position)), nrow(collect.data)))

#############################################
#### Merge with metadata
#############################################

if (is.null(metadata.file) != TRUE){
  meta.sample = read.csv(metadata.file)
}#end if

#lists sample files
sample.names = unique(collect.data$sample)

#combines all the individual samples together
all.samples = c()
for (i in 1:length(sample.names)){

  sample.data = collect.data[collect.data$sample %in% sample.names[i],]

  #Combines metadata if included
  if (is.null(metadata.file) != TRUE){
    sample.data = merge(sample.data, meta.sample, by.x = "sample", by.y = "Sample")
  }#end if

  all.samples = rbind(all.samples, sample.data)
}#end i loop

#############################################
#### Find amino acid changing sites and associate the two
#############################################

sample.names = unique(all.samples$sample)

# Coding sequence of each segment's primary product, built once. Codons are
# indexed into this rather than into the raw segment; for the primary products
# the two happen to coincide, but going through the CDS keeps this script
# honest about where a gene actually stops.
wfabc.cds = list()
for (lc in names(wfabc.reference)){
  for (o in flu_orfs(lc, nchar(wfabc.reference[[lc]])))
    if (o$primary) wfabc.cds[[lc]] = flu_cds_seq(o, wfabc.reference[[lc]])
}

final.data = c()
for (i in 1:length(sample.names)){
  #Subsets to sample data
  sample.data = all.samples[all.samples$sample %in% sample.names[i],]
  #gathers gene names
  gene.names = unique(all.samples$locus)

  #loops through each gene to assess amino acids
  new.gene = c()
  for (j in 1:length(gene.names)){

    gene.data = sample.data[sample.data$locus %in% gene.names[j],]

    if (nrow(gene.data) == 0){ next }

    gene.data$reference_codon = "NA"
    gene.data$alternative_codon = "NA"
    gene.data$reference_aa = "NA"
    gene.data$alternative_aa = "NA"
    gene.data$aa_changing = "NA"

    #Primary-product coding sequence for this segment
    ref.cds = wfabc.cds[[gene.names[j]]]
    if (is.null(ref.cds)){ new.gene = rbind(new.gene, gene.data); next }

    #loops through each row in the gene data to translate
    for (k in 1:nrow(gene.data)){

      #Outside the primary ORF there is no codon to report
      cds.pos = gene.data$cds_position[k]
      aa.pos  = gene.data$aa_position[k]
      if (is.na(cds.pos) || is.na(aa.pos)) next

      #Swap the alternative base in at its CDS index
      alt.cds = ref.cds
      substr(alt.cds, cds.pos, cds.pos) = as.character(gene.data$alternative[k])

      #Subseqs the codons out
      cod.start = (aa.pos - 1) * 3 + 1
      gene.data$reference_codon[k]   = substr(ref.cds, cod.start, cod.start + 2)
      gene.data$alternative_codon[k] = substr(alt.cds, cod.start, cod.start + 2)

      if (nchar(gene.data$reference_codon[k]) < 3) next

      #translates codon
      gene.data$reference_aa[k] = as.character(seqinr::translate(unlist(strsplit(gene.data$reference_codon[k], ""))))
      gene.data$alternative_aa[k] = as.character(seqinr::translate(unlist(strsplit(gene.data$alternative_codon[k], ""))))

      #Checks if it changed the amino acid
      if (gene.data$reference_aa[k] == gene.data$alternative_aa[k]){
        gene.data$aa_changing[k] = "NO"
      } else { gene.data$aa_changing[k] = "YES"}

    }#end k

    new.gene = rbind(new.gene, gene.data)

  }# end j loop

  final.data = rbind(final.data, new.gene)

}#end i

#Save large tab delimited table of all the amino acids
write.table(final.data, paste0(output.directory, "/all_sample_amino_acids.txt"),
            row.names = F, quote = F, sep = "\t")

#############################################
#### Format the allele-frequency trajectories for WFABC
#############################################

# Make sure the metadata columns we need actually exist
needed.cols = c(individual.column, time.column, group.column)
missing.cols = needed.cols[!needed.cols %in% colnames(final.data)]
if (length(missing.cols) != 0){
  stop("Metadata column(s) not found after the merge: ", paste(missing.cols, collapse = ", "),
       ". Check INDIVIDUAL_COLUMN, TIME_COLUMN and GROUP_NAMES in your config file.")
}

#Keep variants that belong to an individual and pass the depth / frequency cutoffs
wfabc.samples = final.data[is.na(final.data[[individual.column]]) != T,]
wfabc.samples = wfabc.samples[wfabc.samples$depth >= min.depth,]
wfabc.samples = wfabc.samples[wfabc.samples$allele_frequency >= min.allele.frequency,]

sample.names = unique(wfabc.samples[[individual.column]])

abc.data = c()
for (i in 1:length(sample.names)){

  sample.data = wfabc.samples[wfabc.samples[[individual.column]] %in% sample.names[i],]

  locus.names = unique(sample.data$locus)

  for (j in 1:length(locus.names)) {

    locus.data = sample.data[sample.data$locus %in% locus.names[j],]

    pos.names = unique(locus.data$position)

    for (k in 1:length(pos.names)){

      # The comma is load-bearing. `df[vec]` on a data.frame selects COLUMNS,
      # not rows, so without it pos.data kept every row and an arbitrary subset
      # of columns: pos.data[[individual.column]] came back NULL and the
      # data.frame below died with "arguments imply differing number of rows:
      # 0, 65". Every other subset in this file has the comma — 286, 355-357,
      # 364, 370 — which is what makes this a typo rather than an intent.
      pos.data = locus.data[locus.data$position %in% pos.names[k], ]

      if (nrow(pos.data) <= 1){
        next
      }

      save.data = data.frame(sample = pos.data[[individual.column]],
                             locus = pos.data$locus,
                             nuc_position = pos.data$position,
                             aa_position = pos.data$aa_position,
                             group = pos.data[[group.column]],
                             time_point = pos.data[[time.column]],
                             sample_size = pos.data$depth,
                             A_allele = round(pos.data$allele_frequency * pos.data$depth, 0))

      if (any(duplicated(save.data$time_point)) == TRUE){
        if (duplicate.handling == "merge"){
          temp.df = aggregate(cbind(sample_size, A_allele) ~ sample+locus+nuc_position+aa_position+group+time_point, data = save.data, sum)
          save.data = temp.df
        }
      }# end if

      if (nrow(save.data) <= 1){
        next
      }

      abc.data = rbind(abc.data, save.data)

    } #end k loop

  }#end j loop

}#end i loop

#############################################
#### Run WFABC and collect the results
#############################################

#Collects the super cool data
header.data = c("sample", "locus", "nuc_position", "aa_position", "group", "number_time_points",
                "Ne_mean", "s_map", "s_lower_bound", "s_upper_bound")

abc.data$generations = abc.data$time_point * generations.per.time

sample.names = unique(abc.data$sample)

#Sets up the results data collection table BEFORE the loop so each value is
#written into the correct column (an earlier version reused the leftover VCF
#variant table here, which shifted every value into the wrong column).
collect.data = data.table::data.table(matrix(as.numeric(0),
                                             nrow = nrow(abc.data),
                                             ncol = length(header.data)))
data.table::setnames(collect.data, header.data)

collect.data[, sample:=as.character(sample)]
collect.data[, locus:=as.character(locus)]
collect.data[, group:=as.character(group)]

# Runs a wfabc command under a timeout. Returns TRUE if it finished, FALSE if it
# was killed for exceeding wfabc.timeout seconds (system() returns 124 / warns).
run_wfabc = function(cmd){
  rc = tryCatch(system(cmd, timeout = wfabc.timeout, ignore.stdout = TRUE, ignore.stderr = TRUE),
                warning = function(w) 124L)
  !identical(as.integer(rc), 124L)
}

# Sites skipped by the near-fixation pre-filter or a wfabc timeout are logged here
skipped.sites = c()

x = 1
for (i in 1:length(sample.names)){

  sample.data = abc.data[abc.data$sample %in% sample.names[i],]

  locus.names = unique(sample.data$locus)

  dir.create(paste0(output.directory, "/", sample.names[i]), showWarnings = FALSE)

  for (j in 1:length(locus.names)) {

    locus.data = sample.data[sample.data$locus %in% locus.names[j],]

    dir.create(paste0(output.directory, "/", sample.names[i], "/", locus.names[j]), showWarnings = FALSE)

    pos.names = unique(locus.data$nuc_position)

    for (k in 1:length(pos.names)){
      print(paste0("sample ", i , " locus ", j, " position ", k))

      pos.data = locus.data[locus.data$nuc_position %in% pos.names[k],]
      pos.data = pos.data[order(pos.data$time_point), ]

      #Clears outliers
      if (nrow(pos.data) >= 3 ){
        pos.data$freq = pos.data$A_allele/pos.data$sample_size
        diff.freq = diff(pos.data$freq)

        # Keep all timepoints by default
        keep <- rep(TRUE, length(pos.data$freq))

        # Flag middle points where the jump is abnormally large in both directions
        # (uses 'm' as the index, not 'x', so the results row counter is not clobbered)
        for (m in 2:(length(pos.data$freq) - 1)) {
          if (abs(pos.data$freq[m] - pos.data$freq[m - 1]) > 0.3 &&
              abs(pos.data$freq[m + 1] - pos.data$freq[m]) > 0.3) {
            keep[m] <- FALSE
          }
        }

        # Always keep the first and last timepoints
        keep[1] <- TRUE
        keep[length(pos.data$freq)] <- TRUE

        pos.data = pos.data[keep,]
      } #end

      #Pre-filter: skip sites where the allele is pinned near one frequency
      #boundary at EVERY time point (no identifiable selection, and the main
      #cause of wfabc_2 hangs). Boundary-crossing trajectories are kept.
      site.tag = paste(sample.names[i], locus.names[j], pos.names[k], sep = "/")
      site.freq = pos.data$A_allele / pos.data$sample_size
      if (all(site.freq >= fixation.cutoff) || all(site.freq <= 1 - fixation.cutoff)){
        print(paste0("  skipping near-fixation site: ", site.tag))
        skipped.sites = rbind(skipped.sites, data.frame(site = site.tag, reason = "near_fixation"))
        next
      }

      dir.create(paste0(output.directory, "/", sample.names[i], "/", locus.names[j], "/", pos.names[k]), showWarnings = FALSE)
      output.file = paste0(output.directory, "/", sample.names[i], "/", locus.names[j], "/", pos.names[k], "/input.txt")

      wfabc_data = c()
      #Header: number of loci, number of time points
      wfabc_data = append(wfabc_data, paste(length(unique(pos.data$locus)),
                                                   length(pos.data$time_point), sep=" "))
      #the time points
      wfabc_data = append(wfabc_data, paste(pos.data$time_point, collapse=","))

      #The total depth and alt allele counts
      total_depths = paste(pos.data$sample_size, collapse = ",")
      alt_allele_counts = paste(pos.data$A_allele, collapse = ",")

      wfabc_data = append(wfabc_data, total_depths)
      wfabc_data = append(wfabc_data, alt_allele_counts)

      # Write the WFABC input file
      writeLines(wfabc_data, output.file)

      #Run WFABC step 1 (estimates Ne)
      setwd(paste0(output.directory, "/", sample.names[i], "/", locus.names[j], "/", pos.names[k]))
      if (!run_wfabc(paste0(wfabc1.bin, " ", output.file))){
        print(paste0("  TIMEOUT in wfabc_1, skipping site: ", site.tag))
        skipped.sites = rbind(skipped.sites, data.frame(site = site.tag, reason = "timeout_wfabc_1"))
        next
      }
      post_N = read.table("input_Ne_bootstrap.txt")

      if (mean(post_N[,1], na.rm = TRUE) <= 0){ next }

      if (is.finite(mean(post_N[,1], na.rm = TRUE))){
        pdf("posterior_density_Ne.pdf", width = 8, height = 8)
        plot(density(post_N[,1]),lwd=2, main= "Ne Posterior", xlab="Ne")
        dev.off()
      } else { next }

      #Run WFABC step 2 (estimates s), widening the prior until it returns output.
      #Clear any stale posterior from a previous run so the first call always runs.
      if (file.exists("input_posterior_s.txt")){ file.remove("input_posterior_s.txt") }
      wfabc2.timed.out = FALSE
      for (bounds in list("-0.3 -max_s 0.3", "-0.3 -max_s 0.3", "-0.1 -max_s 0.1",
                          "-0.01 -max_s 0.01", "-1 -max_s 1")){
        if (file.exists("input_posterior_s.txt") && length(readLines("input_posterior_s.txt")) != 0){ break }
        if (!run_wfabc(paste0(wfabc2.bin, " -ploidy 1 -min_s ", bounds, " ", output.file))){
          wfabc2.timed.out = TRUE; break
        }
      }
      if (wfabc2.timed.out){
        print(paste0("  TIMEOUT in wfabc_2, skipping site: ", site.tag))
        skipped.sites = rbind(skipped.sites, data.frame(site = site.tag, reason = "timeout_wfabc_2"))
        next
      }

      if (!file.exists("input_posterior_s.txt") || length(readLines("input_posterior_s.txt")) == 0){ next }

      #Refine the prior bounds until two consecutive runs agree (KS test, capped)
      ks.p = 0
      counter = 1
      while (ks.p < 0.1 && counter <= max.refine.iter){

        post_s = read.table("input_posterior_s.txt")
        post_s = unlist(as.vector(post_s))

        #estimates lower and upper bound for the next run
        l_bound = quantile(post_s, probs = 0.05)
        u_bound = quantile(post_s, probs = 0.95)

        #Creates a buffer
        range = quantile(post_s, c(0.05, 0.95))
        buffer = 0.05 * diff(range)
        l_bound = l_bound - buffer
        u_bound = u_bound + buffer

        post_s_run1 = post_s

        if (!run_wfabc(paste0(wfabc2.bin, " -ploidy 1 -min_s ", l_bound, " -max_s ", u_bound, " ", output.file))){
          wfabc2.timed.out = TRUE; break
        }

        #This refinement run can occasionally return an empty posterior; if so,
        #stop refining and keep the last valid posterior (post_s_run1).
        if (!file.exists("input_posterior_s.txt") || length(readLines("input_posterior_s.txt")) == 0){
          post_s = post_s_run1
          break
        }

        post_s = read.table("input_posterior_s.txt")
        post_s = unlist(as.vector(post_s))

        post_s_run2 = post_s

        ks.p = suppressWarnings(ks.test(post_s_run1, post_s_run2)$p.value)
        print(paste0("iteration ", counter, " complete."))
        counter = counter+1

      }#end while

      if (wfabc2.timed.out){
        print(paste0("  TIMEOUT in wfabc_2 refinement, skipping site: ", site.tag))
        skipped.sites = rbind(skipped.sites, data.frame(site = site.tag, reason = "timeout_wfabc_2_refine"))
        next
      }

      #Gets the peak (MAP) of the s posterior (post_s is the last valid posterior,
      #held in memory, so this does not re-read a possibly-empty posterior file)
      s_density = density(post_s)
      map_s = s_density$x[which.max(s_density$y)]

      pdf("posterior_density_s.pdf", width = 8, height = 8)
      plot(s_density, main = "Posterior of s with MAP", xlab = "s", lwd=2)
      abline(v = map_s, col = "red", lty = 2)
      text(map_s, max(s_density$y), labels = paste0("MAP = ", round(map_s, 2)), pos = 4, col = "red")
      dev.off()

      s_mcmc = coda::as.mcmc(post_s)
      hpd = coda::HPDinterval(s_mcmc, prob = 0.95)
      lower_bound = hpd[1, "lower"]
      upper_bound = hpd[1, "upper"]

      #Collect the result for this site
      data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = sample.names[i] )
      data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = locus.names[j] )
      data.table::set(collect.data, i = as.integer(x), j = match("nuc_position", header.data), value = pos.names[k] )
      data.table::set(collect.data, i = as.integer(x), j = match("aa_position", header.data), value = unique(pos.data$aa_position) )
      data.table::set(collect.data, i = as.integer(x), j = match("group", header.data), value = unique(pos.data$group) )
      data.table::set(collect.data, i = as.integer(x), j = match("number_time_points", header.data), value = nrow(pos.data) )
      data.table::set(collect.data, i = as.integer(x), j = match("Ne_mean", header.data), value = mean(post_N[,1]) )
      data.table::set(collect.data, i = as.integer(x), j = match("s_map", header.data), value = map_s )
      data.table::set(collect.data, i = as.integer(x), j = match("s_lower_bound", header.data), value = lower_bound )
      data.table::set(collect.data, i = as.integer(x), j = match("s_upper_bound", header.data), value = upper_bound )
      x = x + 1

    }# end k

  }# end j

}# end i

#Record which sites were skipped (near-fixation pre-filter or wfabc timeout)
if (is.null(skipped.sites) != TRUE){
  write.table(skipped.sites, paste0(output.directory, "/wfabc_skipped_sites.txt"),
              row.names = F, quote = F, sep = "\t")
  print(paste0(nrow(skipped.sites), " sites skipped (see wfabc_skipped_sites.txt)"))
}

#Save the WFABC summary table
collect.data = collect.data[collect.data$sample != 0,]
save.data = collect.data

write.table(save.data, paste0(output.directory, "/wfabc_summary.txt"),
            row.names = F, quote = F, sep = "\t")

#############################################
#### Frequency Increment Test (drift vs selection)
#############################################
# WFABC answers "what selection coefficient best explains this trajectory".
# The FIT (Feder, Kryazhimskiy & Plotkin 2014) answers the prior question:
# is the trajectory distinguishable from DRIFT at all? The two are reported
# side by side because they can disagree, and when they do the disagreement is
# the finding — a large s_map whose credible interval spans zero and whose FIT
# p-value is unremarkable is drift wearing a selection coefficient.
#
# Ported from analyze_drift.R (Cow_Comparison), generalised: no hardcoded paths
# and the group comes from GROUP_NAMES rather than a study-specific vaccine
# label. Column names are kept identical to that script's FIT_results.csv so
# the two analyses stay directly comparable.
#
# It reads the per-site input.txt files WFABC already wrote — line 2 is the time
# points, line 3 the depths, line 4 the alt counts.
fit.one = function(f, root){
  L = readLines(f, warn = FALSE)
  if (length(L) < 4) return(NULL)
  tp  = as.numeric(strsplit(L[2], ",")[[1]])
  dep = as.numeric(strsplit(L[3], ",")[[1]])
  alt = as.numeric(strsplit(L[4], ",")[[1]])
  # FIT needs at least three time points: the test is a t-test over the
  # INCREMENTS, so two time points give a single increment and no variance.
  if (length(tp) < 3) return(NULL)
  # (alt+0.5)/(dep+1) keeps the frequency strictly inside (0,1); at 0 or 1 the
  # rescaling below divides by zero.
  v = (alt + 0.5) / (dep + 1)
  o = order(tp); tp = tp[o]; v = v[o]
  dt = diff(tp); dv = diff(v)
  Y = dv / sqrt(2 * head(v, -1) * (1 - head(v, -1)) * dt)
  Y = Y[is.finite(Y)]
  if (length(Y) < 2 || sd(Y) == 0) return(NULL)
  tt = t.test(Y)
  parts = strsplit(sub(paste0(root, "/"), "", f), "/")[[1]]   # animal/locus/pos
  if (length(parts) < 3) return(NULL)
  data.frame(animal = parts[1], locus = parts[2], pos = as.numeric(parts[3]),
             n_tp = length(tp), FIT_t = unname(tt$statistic), FIT_p = tt$p.value,
             mean_incr = mean(Y), stringsAsFactors = FALSE)
}

traj.files = list.files(output.directory, pattern = "^input.txt$",
                        recursive = TRUE, full.names = TRUE)
fit = NULL
if (length(traj.files)){
  fit = do.call(rbind, lapply(traj.files, fit.one, root = output.directory))
}
if (is.null(fit) || nrow(fit) == 0){
  print(paste0("FIT: no trajectory had 3+ time points (", length(traj.files),
               " trajectories) — FIT_results.csv not written"))
} else {
  fit$FIT_fdr = p.adjust(fit$FIT_p, "BH")

  # Join WFABC's own estimates on. wfabc_dir is a LOGICAL, not a path: TRUE when
  # the ABC credible interval excludes zero. The name is inherited from
  # analyze_drift.R and kept for compatibility even though it reads like a
  # directory.
  # as.data.frame is not decoration: collect.data is a data.table, where
  # sm[, keep, drop=FALSE] treats `keep` as a column NAME rather than a vector
  # of names and dies with "column name 'keep' is not found". Converting once
  # here means the rest of this block is plain data.frame semantics regardless
  # of what the pipeline upstream happens to hand over.
  sm = as.data.frame(save.data)
  sm$key  = paste(sm$sample, sm$locus, sm$nuc_position)
  fit$key = paste(fit$animal, fit$locus, fit$pos)
  keep = intersect(c("key","s_map","s_lower_bound","s_upper_bound","Ne_mean"), names(sm))
  fit = merge(fit, sm[, keep, drop = FALSE], by = "key", all.x = TRUE)
  fit$wfabc_dir = !is.na(fit$s_lower_bound) &
                  (fit$s_lower_bound > 0 | fit$s_upper_bound < 0)

  # Group label, when the metadata carried one.
  if (group.column %in% names(all.samples)){
    grp.map = as.data.frame(all.samples)          # same data.table caution
    grp.map = grp.map[!duplicated(grp.map[[individual.column]]), ]
    fit$grp = grp.map[[group.column]][match(fit$animal,
                                            as.character(grp.map[[individual.column]]))]
  } else {
    fit$grp = NA
  }

  ord = c("grp","animal","locus","pos","n_tp","FIT_t","FIT_p","FIT_fdr",
          "mean_incr","s_map","wfabc_dir","Ne_mean")
  fit = fit[, intersect(ord, names(fit)), drop = FALSE]
  write.csv(fit, paste0(output.directory, "/FIT_results.csv"), row.names = FALSE)

  print(paste0("FIT: ", nrow(fit), " of ", length(traj.files),
               " trajectories had 3+ time points; ",
               sum(fit$FIT_fdr < 0.05, na.rm = TRUE), " significant at BH-FDR<0.05; ",
               sum(fit$n_tp >= 4, na.rm = TRUE), " with 4+ time points (df>=2)"))
}

#############################################
#### Summary figures (generic + batch-safe)
#############################################
# Every figure is written to a PDF in output.directory. The whole section is
# guarded so the script never fails at the end: if ggplot2 is missing, a table
# is empty, or a single plot errors, it prints a message and moves on. It uses
# only the normalized columns the pipeline always produces (save.data for the
# selection results, abc.data for the allele-frequency trajectories), so it
# works for any dataset without dataset-specific column names.

if (!requireNamespace("ggplot2", quietly = TRUE)){
  message("ggplot2 not installed; skipping summary figures.")
} else {
  library(ggplot2)

  sel  = as.data.frame(save.data)
  traj = if (is.null(abc.data)) NULL else as.data.frame(abc.data)
  if (!is.null(traj) && nrow(traj) > 0){ traj$frequency = traj$A_allele / traj$sample_size }

  # Significance: a site is flagged when its 95% HPD interval for s excludes 0
  # (credible interval entirely above 0 = positive, or entirely below 0 = negative).
  # This is a per-site Bayesian criterion with no multiple-testing correction.
  sel$sig = !is.na(sel$s_lower_bound) & !is.na(sel$s_upper_bound) &
            (sel$s_lower_bound > 0 | sel$s_upper_bound < 0)

  # Colour-blind-friendly default palette (Okabe-Ito), keyed to whatever groups
  # are present so all figures share one consistent scheme for any dataset. The
  # palette recycles if a dataset has more than eight groups.
  group.palette.base = c("#E69F00", "#CC79A7", "#0072B2", "#009E73",
                         "#D55E00", "#56B4E9", "#F0E442", "#000000")
  all.groups = sort(unique(c(as.character(sel$group),
                             if (!is.null(traj)) as.character(traj$group))))
  group.cols = stats::setNames(rep(group.palette.base, length.out = length(all.groups)),
                               all.groups)

  save_plot = function(p, file, w = 12, h = 8){
    tryCatch(ggsave(file.path(output.directory, file), p, width = w, height = h),
             error = function(e) message("could not save ", file, ": ", conditionMessage(e)))
  }

  ## 1. Selection coefficient (s_map) with 95% HPD, per locus and position, by group
  tryCatch({
    if (nrow(sel) > 0){
      p1 = ggplot(sel, aes(x = factor(aa_position), y = s_map, colour = factor(group))) +
        geom_hline(yintercept = 0, colour = "black") +
        geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted", colour = "grey50") +
        geom_errorbar(aes(ymin = s_lower_bound, ymax = s_upper_bound),
                      width = 0.25, position = position_dodge(width = 0.6)) +
        geom_point(position = position_dodge(width = 0.6), size = 2) +
        # asterisk above sites whose 95% HPD interval excludes 0 (significant)
        geom_text(data = sel[sel$sig, , drop = FALSE],
                  aes(y = s_upper_bound + 0.03, label = "*", group = factor(group)),
                  position = position_dodge(width = 0.6),
                  colour = "black", size = 5, fontface = "bold", show.legend = FALSE) +
        facet_wrap(~ locus, scales = "free_x") +
        scale_colour_manual(values = group.cols, drop = FALSE) +
        labs(x = "Amino acid position", y = "s (MAP +/- 95% HPD)", colour = "Group",
             title = "Selection coefficients per locus and position",
             subtitle = "* = 95% HPD interval excludes 0 (uncorrected)") +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(p1, "plot_selection_coefficients.pdf")
    }
  }, error = function(e) message("plot 1 (selection coefficients) failed: ", conditionMessage(e)))

  ## 2. Number of sites under strong selection (|s| >= 0.1) per locus and group
  tryCatch({
    if (nrow(sel) > 0){
      sel$direction = ifelse(sel$s_map >=  0.1, "positive (s >= 0.1)",
                      ifelse(sel$s_map <= -0.1, "negative (s <= -0.1)", "neutral (|s| < 0.1)"))
      cnt = as.data.frame(table(locus = sel$locus, group = sel$group, direction = sel$direction))
      cnt = cnt[cnt$Freq > 0, ]
      write.table(cnt, file.path(output.directory, "selection_counts.txt"),
                  row.names = FALSE, quote = FALSE, sep = "\t")
      strong = cnt[cnt$direction != "neutral (|s| < 0.1)", ]
      if (nrow(strong) > 0){
        p2 = ggplot(strong, aes(x = locus, y = Freq, fill = factor(group))) +
          geom_col(position = position_dodge(width = 0.8)) +
          facet_wrap(~ direction) +
          scale_fill_manual(values = group.cols, drop = FALSE) +
          labs(x = "Locus", y = "Number of sites", fill = "Group",
               title = "Sites under strong selection (|s| >= 0.1)") +
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_plot(p2, "plot_selection_counts.pdf")
      }
    }
  }, error = function(e) message("plot 2 (selection counts) failed: ", conditionMessage(e)))

  ## 3. Allele-frequency trajectories through time, per locus
  tryCatch({
    if (!is.null(traj) && nrow(traj) > 0){
      traj$trajectory = paste(traj$sample, traj$nuc_position, sep = "_")
      p3 = ggplot(traj, aes(x = time_point, y = frequency,
                            group = trajectory, colour = factor(group))) +
        geom_line(alpha = 0.5) + geom_point(size = 1) +
        facet_wrap(~ locus, scales = "free_x") +
        scale_colour_manual(values = group.cols, drop = FALSE) +
        labs(x = "Time point", y = "Alternative allele frequency", colour = "Group",
             title = "Allele-frequency trajectories through time") +
        theme_bw()
      save_plot(p3, "plot_frequency_trajectories.pdf")
    }
  }, error = function(e) message("plot 3 (frequency trajectories) failed: ", conditionMessage(e)))

  ## 4. Allele frequency by amino-acid position across each locus
  tryCatch({
    if (!is.null(traj) && nrow(traj) > 0){
      p4 = ggplot(traj, aes(x = aa_position, y = frequency, colour = factor(group))) +
        geom_point(alpha = 0.6, size = 1.5) +
        facet_wrap(~ locus, scales = "free_x") +
        scale_colour_manual(values = group.cols, drop = FALSE) +
        labs(x = "Amino acid position", y = "Alternative allele frequency", colour = "Group",
             title = "Allele frequency by position") +
        theme_bw()
      save_plot(p4, "plot_allele_frequency_by_position.pdf")
    }
  }, error = function(e) message("plot 4 (allele frequency by position) failed: ", conditionMessage(e)))
}

#########################
###### END SCRIPT
#########################
