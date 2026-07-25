#### Prior to this script, run the Flumina pipeline to obtain variant calls in VCF files
#### Next, run 1_convertVCFtoTable.R to get a vcf table, and place output in "vcftable.path" below 

# Required packages:
# 1. data.table
# 2. foreach
# 3. doParallel
# 4. snow
# 5. Biostrings
# 6. seqinr

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

#Get multifile databases together
reference = Biostrings::readDNAStringSet(reference.path, format = "fasta")
vcf.data = read.csv(vcftable.path, header = T)
sample.names = unique(vcf.data$sample)

# Sets up multiprocessing
my.cluster = parallel::makeCluster(threads, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

#Loop for cd-hit est reductions
foreach::foreach(i = seq_along(sample.names), .packages = c("foreach", "Biostrings", "seqinr", "data.table")) %dopar% {
  
#for (i in 1:length(sample.names)){
  #Subsets to sample data
  sample.data = vcf.data[vcf.data$sample %in% sample.names[i],]
  #gathers gene names
  gene.names = unique(vcf.data$locus)
  
  #loops through each gene to assess amino acids
  new.gene = c()
  for (j in 1:length(gene.names)){
    
    #Creates empty spots for new variables
    gene.data = sample.data[sample.data$locus %in% gene.names[j],]
    
    if (nrow(gene.data) == 0){ next }
    
    gene.data$reference_codon = "NA"
    gene.data$alternative_codon = "NA"
    gene.data$reference_aa = "NA"
    gene.data$alternative_aa = "NA"
    gene.data$aa_changing = "NA"
    
    #Subsets to reference for specific gene, translates
    ref.seq = as.character(reference[names(reference) == gene.names[j]])

    #loops through each row in the gene data to translate
    for (k in 1:nrow(gene.data)){
    
      #skips if no data for gene  
      if (nrow(gene.data) == 0){ next }
      
      #Converts gene reference to character vector to change to alternative allele position
      gene.char = unlist(strsplit(as.character(ref.seq), ""))
      gene.char[gene.data$position[k]] = gene.data$alternative[k]
      new.seq = Biostrings::DNAStringSet(paste(gene.char, collapse = ""))

      #Subseqs the codons out
      gene.data$reference_codon[k] = Biostrings::subseq(ref.seq, 
                                                        start = (gene.data$aa_position[k] - 1) * 3 + 1, 
                                                        end = (gene.data$aa_position[k] - 1) * 3 + 3 ) 
      
      gene.data$alternative_codon[k] = Biostrings::subseq(new.seq, 
                                                          start = (gene.data$aa_position[k] - 1) * 3 + 1, 
                                                          end = (gene.data$aa_position[k] - 1) * 3 + 3 ) 
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
for (i in 1:length(aa.sample)){
  
  sample.data = read.csv(paste0(output.directory, "/aa_db/", aa.sample[i]), header = TRUE, sep = ",")
  
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
    dup.key = do.call(paste, c(all.samples[, c(keys, "locus", "position", "alternative")],
                               sep = "|"))
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




