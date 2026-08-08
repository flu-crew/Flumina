#### Prior to this script, run the Flumina pipeline to obtain variant calls in VCF files
#### Next, run 1_convertVCFtoTable.R to get a vcf table, and place output in "vcftable.path" below 

# Required packages:
# 1. data.table
# 2. foreach
# 3. doParallel
# 4. snow
# 5. seqinr

#Takes 5 minutes to run on 500 samples

args = commandArgs(trailingOnly = TRUE)
#args = "config.cfg"

#args = "config.cfg"
# Function to read and parse configuration file
lines <- readLines(args)

#makes a list and loads stuff in with an equal sign
config <- list()
for (line in lines) {
  line <- trimws(line)  # Remove leading and trailing whitespaces
  if (nchar(line) != 0 && !startsWith(line, "#")) {  # Ignore empty lines and comments
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      value <- trimws(parts[2])
      config[[key]] <- value
    }# end if
  }#end if
}#end for

#Define these

# output.directory used from 1_convertVCFtoTable.R
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/variant_analysis")

# file name (or path if moved) of the table from 1_convertVCFtoTable.R
vcftable.path = paste0(output.directory, "/variant-table.csv")

#Single reference used for variant calling, full path if not in working directory
reference.path = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY),"/reference.fa")

#Optional for merging metadata with AA data, set to NULL if none available
metadata.file = gsub("\"", "", config$METADATA)
if(length(metadata.file) == 0L || metadata.file == "NULL") {
  metadata.file <- NULL
}

#Set multithreading and memory usage
threads = as.numeric(gsub("\"", "", config$THREADS))

# Config helper with defaults (returns `default` when the key is missing/blank/NULL).
# Robust to interactive runs where the config file was not parsed into `config`.
if (!exists("config") || !is.list(config)) config = list()
cfg = function(key, default = NULL){
  v = config[[key]]
  if (is.null(v)) return(default)
  v = gsub("\"", "", trimws(v))
  if (!nzchar(v) || v == "NULL") return(default)
  v
}

# Depth-artifact guards applied to the combined amino-acid table (see below).
# MIN_DEPTH defaults to 100x: at low template input LoFreq/GATK4 report false
# fixations (founder/jackpot effect), so a real depth floor is enforced by
# default even if a study omits the key. MIN_QUALITY defaults off (0).
# DEDUP_KEYS is opt-in: a comma-separated list of metadata columns identifying
# one biological sample (e.g. "Animal.ID,DPI,Quarter"); when set, replicate
# libraries of the same sample are collapsed by keeping the deepest call per
# sample x locus x position x alternative. Unset -> no dedup (default).
min.depth   = as.numeric(cfg("MIN_DEPTH", "100"))
min.quality = as.numeric(cfg("MIN_QUALITY", "0"))
dedup.keys  = cfg("DEDUP_KEYS", "")

# output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Bird_Flu/bird_flu_new/variant_analysis"
# vcftable.path = paste0(output.directory, "/variant-table.csv")
# reference.path = paste0("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Bird_Flu/bird_flu_new/Reference/reference.fa")
# threads = 4
# metadata.file = NULL

#############################################
#### Should not need to modify below here
#############################################

#Sets working directory and creates output
dir.create(paste0(output.directory, "/aa_db"))
require(foreach)

# Coding intervals for every product, shared with convertVCFtoTable.R so the
# two cannot disagree about where a gene starts.
script.dir = dirname(sub("^--file=", "",
                         grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(script.dir) || !nzchar(script.dir)) script.dir = "."
source(file.path(script.dir, "fluORFs.R"))

#Get multifile databases together
reference = flu_read_fasta(reference.path)

# na.strings = "" is NOT optional here. The neuraminidase gene is called "NA",
# so read.csv's default turns every neuraminidase row's `product` (and, once
# the locus is shortened, its `locus`) into a missing value and the entire
# segment silently disappears — 797 of 18,771 rows on the swine WGS data.
# outputSummary.R reads the amino-acid table with the same setting for exactly
# this reason.
vcf.data = read.csv(vcftable.path, header = TRUE, na.strings = "")
sample.names = unique(vcf.data$sample)

# Older variant tables have no `product` column. Rather than silently translate
# them in the wrong frame, annotate them here the same way convertVCFtoTable.R
# now does.
if (!"product" %in% colnames(vcf.data)){
  message("variant-table.csv predates per-product annotation; annotating now.")
  vcf.data = flu_annotate_positions(vcf.data, reference,
                                    locus.col = "locus", pos.col = "position")
}

# The SPLICED coding sequence of every product, built once.
#
# This is the correction that matters. The old code took codons straight out of
# the segment at (aa_position-1)*3+1, which silently assumes the CDS is
# contiguous and starts at nucleotide 1. For M2, NEP and PA-X it is neither, so
# every codon it produced for those products was read out of the wrong place —
# and past the end of M1 and NS1 it happily translated through the stop codon
# and returned a plausible-looking amino acid for a residue that does not exist.
# Indexing into the spliced CDS instead makes the junction a non-event.
cds.lookup = list()
for (lc in names(reference)){
  for (o in flu_orfs(lc, nchar(reference[[lc]])))
    cds.lookup[[paste(lc, o$gene, sep = "|")]] = flu_cds_seq(o, reference[[lc]])
}

# One translation unit per (segment, product) pair present in the data
gene.keys = unique(vcf.data[, c("locus", "product")])
gene.keys = gene.keys[!is.na(gene.keys$product), , drop = FALSE]

# Sets up multiprocessing
my.cluster = parallel::makeCluster(threads, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

#Loop for cd-hit est reductions
foreach::foreach(i = seq_along(sample.names), .packages = c("foreach", "seqinr", "data.table")) %dopar% {

#for (i in 1:length(sample.names)){
  #Subsets to sample data
  sample.data = vcf.data[vcf.data$sample %in% sample.names[i],]

  #loops through each (segment, product) pair to assess amino acids
  new.gene = c()
  for (j in seq_len(nrow(gene.keys))){

    #Creates empty spots for new variables
    gene.data = sample.data[sample.data$locus   %in% gene.keys$locus[j] &
                            sample.data$product %in% gene.keys$product[j], ]

    if (nrow(gene.data) == 0){ next }

    gene.data$reference_codon = "NA"
    gene.data$alternative_codon = "NA"
    gene.data$reference_aa = "NA"
    gene.data$alternative_aa = "NA"
    gene.data$aa_changing = "NA"

    #The spliced CDS of this product; codons are indexed into THIS, not the segment
    ref.cds = cds.lookup[[paste(gene.keys$locus[j], gene.keys$product[j], sep = "|")]]
    if (is.null(ref.cds)){ new.gene = rbind(new.gene, gene.data); next }

    #loops through each row in the gene data to translate
    for (k in 1:nrow(gene.data)){

      #A position in no coding region has no codon. Leave it "NA" rather than
      #inventing one — this is exactly the case the old code got wrong.
      cds.pos = gene.data$cds_position[k]
      aa.pos  = gene.data$aa_position[k]
      if (is.na(cds.pos) || is.na(aa.pos)) next

      #Swap the alternative base in at its CDS index
      alt.cds = ref.cds
      substr(alt.cds, cds.pos, cds.pos) = as.character(gene.data$alternative[k])

      #Subseqs the codons out of the spliced CDS
      cod.start = (aa.pos - 1) * 3 + 1
      gene.data$reference_codon[k]   = substr(ref.cds, cod.start, cod.start + 2)
      gene.data$alternative_codon[k] = substr(alt.cds, cod.start, cod.start + 2)

      #A truncated final codon cannot be translated
      if (nchar(gene.data$reference_codon[k]) < 3) next

      #translates codon
      gene.data$reference_aa[k] = as.character(seqinr::translate(unlist(strsplit(gene.data$reference_codon[k], ""))))
      gene.data$alternative_aa[k] = as.character(seqinr::translate(unlist(strsplit(gene.data$alternative_codon[k], ""))))

      #Checks if it changed the amino acid
      if (gene.data$reference_aa[k] == gene.data$alternative_aa[k]){
        gene.data$aa_changing[k] = "NO"
      } else { gene.data$aa_changing[k] = "YES"}

    }#end k

    #Saves all data
    new.gene = rbind(new.gene, gene.data)

  }# end j loop

  #Writes data for each sample
  write.csv(new.gene, paste0(output.directory, "/aa_db/", sample.names[i], ".csv"),
            row.names = F, quote = F)

}#end i
  
parallel::stopCluster(cl = my.cluster)


#########################
###### Make a combined spreadsheet
#########################

# if there is metadata to join
if (is.null(metadata.file) != TRUE){
  meta.sample = read.csv(metadata.file)
}#end if

#lists sample files
aa.sample = list.files(paste0(output.directory, "/aa_db"))

#combines all the individual samples together
all.samples = c()
for (i in seq_along(aa.sample)){
  
  # na.strings = "" again — see the note above; "NA" is a gene name here
  sample.data = read.csv(paste0(output.directory, "/aa_db/", aa.sample[i]),
                         header = TRUE, sep = ",", na.strings = "")
  
  #Combines metadata if included
  if (is.null(metadata.file) != TRUE){
    sample.data = merge(sample.data, meta.sample, by.x = "sample", by.y = "Sample")
  }#end if
  
  sample.data$locus = gsub("^A_", "", sample.data$locus)
  sample.data$locus = gsub("_[A-Z][0-9]+$", "", sample.data$locus)
  all.samples = rbind(all.samples, sample.data)
}#end i loop

#############################################
#### Depth-artifact guards (depth/quality floor + optional swab dedup)
#############################################
# Low-input libraries produce spurious apparent fixations (founder/jackpot
# effect); enforce a real read-depth floor here so the amino-acid table (and
# every downstream high-frequency/curated analysis that reads it) is clean.
all.samples$depth   = suppressWarnings(as.numeric(all.samples$depth))
all.samples$quality = suppressWarnings(as.numeric(all.samples$quality))
n0 = nrow(all.samples)
all.samples = all.samples[!is.na(all.samples$depth) & all.samples$depth >= min.depth, ]
if (min.quality > 0){
  all.samples = all.samples[!is.na(all.samples$quality) & all.samples$quality >= min.quality, ]
}
cat(sprintf("Depth/quality filter (MIN_DEPTH=%s, MIN_QUALITY=%s): %d -> %d calls\n",
            min.depth, min.quality, n0, nrow(all.samples)))

# Optional swab deduplication: collapse replicate libraries of the same
# biological sample by keeping the DEEPEST call per DEDUP_KEYS x locus x
# position x alternative. Disabled unless DEDUP_KEYS is set in the config.
if (nzchar(dedup.keys)){
  keys = trimws(strsplit(dedup.keys, ",")[[1]])
  miss = setdiff(keys, colnames(all.samples))
  if (length(miss) > 0){
    warning("DEDUP_KEYS column(s) not found (", paste(miss, collapse = ", "),
            "); skipping swab dedup.")
  } else {
    n1 = nrow(all.samples)
    all.samples = all.samples[order(-all.samples$depth), ]
    # `product` belongs in the key: one nucleotide can code in two reading
    # frames, so (locus, position, alternative) alone would collapse a PA row
    # and its PA-X counterpart into one and silently drop a real annotation.
    dedup.cols = c(keys, "locus", "product", "position", "alternative")
    dedup.cols = dedup.cols[dedup.cols %in% colnames(all.samples)]
    dup.key = do.call(paste, c(all.samples[, dedup.cols], sep = "|"))
    all.samples = all.samples[!duplicated(dup.key), ]
    cat(sprintf("Swab dedup on [%s + locus/position/alternative]: %d -> %d calls\n",
                paste(keys, collapse = ","), n1, nrow(all.samples)))
  }
}

#Save large tab delimited table of all the amino acids
write.table(all.samples, paste0(output.directory, "/all_sample_amino_acids.txt"),
            row.names = F, quote = F, sep = "\t")


#########################
###### end script
#########################




