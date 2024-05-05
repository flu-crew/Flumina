
#Debegging

args = commandArgs(trailingOnly = TRUE)

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
overwrite = is.logical(gsub("\"", "", config$OVERWRITE))
irma.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/IRMA_results")
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/IRMA-consensus-contigs")


#Quick checks
if (is.null(irma.directory) == TRUE){ stop("Please provide the read directory.") }
if (file.exists(irma.directory) == F){ stop("read folder not found.") }

# Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
if (dir.exists(output.directory) == F) {
  dir.create(output.directory)
} else {
  if (overwrite == TRUE) {
    system(paste0("rm -r ", output.directory))
    dir.create(output.directory)
  }
} # end else

#Read in sample data
sample.names = list.files(irma.directory, recursive = F, full.names = F)

#Loops through each sample and  moves them to new directory
for (i in 1:length(sample.names)) {
  #################################################
  ### Part A: prepare for loading and checks
  #################################################
  #Gets the reads for the sample
  sample.files = list.files(paste0(irma.directory, "/", sample.names[i]), full.names = T)
  fasta.files = sample.files[grep(".fasta$", sample.files)]
  
  #check
  if (file.exists(paste0(output.directory, "/", sample.names[i], ".fasta")) == TRUE){
    unlink(paste0(output.directory, "/", sample.names[i], ".fasta"))
  }
  
  #Starts to finally look for Haplotypes! *here
  system(paste0(
    "cat ", paste0(fasta.files, collapse = " "),
    " > ", output.directory, "/", sample.names[i], ".fasta"
  ))



}# end i loop


# END SCRIPT
