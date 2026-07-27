#### Prior to this script, run the Flumina pipeline to obtain variant calls in VCF files

# Required packages:
# 1. data.table package

#Takes 2-3 minutes to run on 500 samples

args = commandArgs(trailingOnly = TRUE)
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

#output directory for analysis
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/variant_analysis")

#vcf directory name, full path if not in working directory
vcf.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/vcf_files")

#name for the table
save.name = "variant-table"

# Amino-acid positions come from the actual coding intervals, not from
# ceiling(position/3) — see Scripts/fluORFs.R for why that is wrong for M2,
# NEP, PA-X and PB1-F2. The ORF definitions are derived from the reference
# FASTA, NOT from reference_gtf/, because that directory is written later by
# the optional SNPGenie step and does not exist when this script runs.
script.dir = dirname(sub("^--file=", "",
                         grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(script.dir) || !nzchar(script.dir)) script.dir = "."
source(file.path(script.dir, "fluORFs.R"))

reference.path = gsub("\"", "", config$REFERENCE_FILE)
if (length(reference.path) == 0L || !nzchar(reference.path) ||
    reference.path == "NULL" || !file.exists(reference.path)) {
  reference.path = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/reference.fa")
}

#output.directory = "/Volumes/Extreme_SSD/Bailey_project/variant_analysis"
#vcf.directory = "/Volumes/Extreme_SSD/Bailey_project/vcf_files"


#############################################
#### Should not need to modify below here
#############################################

#### LoFreq

#the string or name of the VCF file for data anaylsis 
vcf.string = "lofreq-called-variants.vcf" #or "gatk4-filtered-snps.vcf"

#Creates output directory
dir.create(output.directory)

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
    #Collect data
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    #Sample data
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = 0 )
    #Length data
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
    
    #Collect data
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    #Sample data
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = VCF$'#CHROM'[j] )
    #Length data
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = VCF$POS[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = VCF$REF[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = VCF$ALT[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = VCF$QUAL[j] )
    
    #find depth
    depth.val = as.numeric(gsub(";", "", gsub(";.*", "", gsub(".*DP=", "", VCF[j,]$INFO))) ) 
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = depth.val)
    
    #find depth
    if (length(grep("MQ=", VCF[j,]$INFO)) != 0){
      mq.val = as.numeric(gsub(";.*", "", gsub(".*;MQ=", "", VCF[j,]$INFO)))
    } else { mq.val = NA }
    
    #Map quality
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = mq.val)

    #find depth
    freq.val = as.numeric(gsub(";.*", "", gsub(".*;AF=", "", VCF[j,]$INFO)))
    data.table::set(collect.data, i = as.integer(x), j = match("allele_frequency", header.data), value = freq.val)
    
    # aa_position is NOT computed here any more. It depends on which product a
    # position codes for, and a position can code for two, so it is filled in
    # after both callers have been read — see the annotation block below.
    #counter goes counting
    x = x + 1
  }#end j loop
  
}#end i loop

#Removes empty samples
collect.data = collect.data[collect.data$sample != 0,]
collect.data = collect.data[collect.data$locus != 0,]

#Sometimes the amino acid T will turn into TRUE
collect.data$alternative[collect.data$alternative == "TRUE"] = "T"
collect.data$reference[collect.data$reference == "TRUE"] = "T"

lofreq.data = collect.data

#############################################
#### GATK4
#############################################

#the string or name of the VCF file for data anaylsis 
vcf.string = "gatk4-filtered-snps.vcf" #or "gatk4-filtered-snps.vcf"

#Creates output directory
dir.create(output.directory)

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
    
    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "GATK4")
    #Collect data
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    #Sample data
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = 0 )
    #Length data
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
    
    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "GATK4")
    
    #Collect data
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    #Sample data
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = VCF$'#CHROM'[j] )
    #Length data
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = VCF$POS[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = VCF$REF[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = VCF$ALT[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = VCF$QUAL[j] )
    
    #find depth
    depth.val = as.numeric(gsub(";", "", gsub(";.*", "", gsub(".*DP=", "", VCF[j,]$INFO))) ) 
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = depth.val)
    
    #find depth
    if (length(grep("MQ=", VCF[j,]$INFO)) != 0){
      mq.val = as.numeric(gsub(";.*", "", gsub(".*;MQ=", "", VCF[j,]$INFO)))
    } else { mq.val = NA }
    
    #Map quality
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = mq.val)
    
    #find depth
    freq.val = as.numeric(gsub(";.*", "", gsub(".*;AF=", "", VCF[j,]$INFO)))
    data.table::set(collect.data, i = as.integer(x), j = match("allele_frequency", header.data), value = freq.val)
    
    # aa_position is NOT computed here any more. It depends on which product a
    # position codes for, and a position can code for two, so it is filled in
    # after both callers have been read — see the annotation block below.
    #counter goes counting
    x = x + 1
  }#end j loop
  
}#end i loop

#Removes empty samples
collect.data = collect.data[collect.data$sample != 0,]
collect.data = collect.data[collect.data$locus != 0,]

#Sometimes the amino acid T will turn into TRUE
collect.data$alternative[collect.data$alternative == "TRUE"] = "T"
collect.data$reference[collect.data$reference == "TRUE"] = "T"

final.data = rbind(lofreq.data, collect.data)

#############################################
#### Amino-acid annotation, through the real coding intervals
#############################################
# Eight segments, twelve proteins. A nucleotide inside PA may also code for
# PA-X in a different frame; one inside MP past nucleotide 756 codes for M2 and
# for nothing else. So this expands to ONE ROW PER PRODUCT rather than assuming
# a single reading frame starting at nucleotide 1.
#
# `locus` deliberately stays the SEGMENT. Everything downstream that keys on it
# — the curated-database join in outputSummary.R above all — keeps working
# unchanged; use the new `product` column to separate reading frames.

final.data = as.data.frame(final.data, stringsAsFactors = FALSE)
n.before   = nrow(final.data)

reference = flu_read_fasta(reference.path)
cat("Annotating amino-acid positions from", reference.path, "\n")
print(flu_orf_table(names(reference), nchar(reference)))

final.data = flu_annotate_positions(final.data, reference,
                                    locus.col = "locus", pos.col = "position")

prim = sum(final.data$product_primary, na.rm = TRUE)
cat(sprintf("Amino-acid annotation: %d calls -> %d rows (%d primary, %d secondary ORF)\n",
            n.before, nrow(final.data), prim, nrow(final.data) - prim))
if (nrow(final.data) < n.before)
  cat(sprintf("  %d call(s) fell in no coding region (UTR, stop codon or intron) and were dropped\n",
              n.before - nrow(final.data)))
tab = table(final.data$product)
cat("  rows per product:", paste(names(tab), tab, sep = "=", collapse = "  "), "\n")

#Saves the data
write.csv(final.data, paste0(output.directory, "/", save.name, ".csv"),
          row.names = F, quote = F)




#########################
###### END SCRIPT
#########################














