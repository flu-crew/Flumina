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
    
    #Finds amino acid positions from nucleotide positions
    aa.pos = ceiling(VCF$POS[j]/3)
    data.table::set(collect.data, i = as.integer(x), j = match("aa_position", header.data), value = aa.pos)
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
    
    #Finds amino acid positions from nucleotide positions
    aa.pos = ceiling(VCF$POS[j]/3)
    data.table::set(collect.data, i = as.integer(x), j = match("aa_position", header.data), value = aa.pos)
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

#Saves the data
write.csv(final.data, paste0(output.directory, "/", save.name, ".csv"),
          row.names = F, quote = F)




#########################
###### END SCRIPT
#########################














