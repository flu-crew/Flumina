
# Debugging.

args = commandArgs(trailingOnly = TRUE)

# Function to read and parse configuration file
lines <- readLines(args)

# Parse key-value pairs into a named list.
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
  } # End of if block.
} # End of for loop.

# Define these parameters.
# `is.logical()` on a character vector is always FALSE, so OVERWRITE would not
# take effect. Parse the string as a logical value instead.
overwrite.raw = gsub("\"", "", config$OVERWRITE)
overwrite = length(overwrite.raw) > 0L &&
            toupper(trimws(overwrite.raw)) %in% c("TRUE", "T", "YES", "1")
irma.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/IRMA_results")
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/IRMA-consensus-contigs")
# IRMA amended consensus (N-masked, IUPAC ambiguity codes), one FASTA per sample.
amended.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/IRMA-amended-contigs")


# Perform basic validation.
if (is.null(irma.directory) == TRUE){ stop("Please provide the read directory.") }
if (file.exists(irma.directory) == F){ stop("read folder not found.") }

# Create the output directories. Overwrite them when requested.
for (out.dir in c(output.directory, amended.directory)) {
  if (dir.exists(out.dir) == F) {
    dir.create(out.dir)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", out.dir))
      dir.create(out.dir)
    }
  } # end else
} # end out.dir loop

# Read the sample data.
sample.names = list.files(irma.directory, recursive = F, full.names = F)

# Iterate over samples and move their files to the destination directory.
for (i in seq_along(sample.names)) {
  #################################################
  ### Part A: prepare inputs and perform checks
  #################################################
  # Locate the reads for the sample.
  sample.files = list.files(paste0(irma.directory, "/", sample.names[i]), full.names = T)
  fasta.files = sample.files[grep(".fasta$", sample.files)]
  
  # Remove the destination directory when it is empty.
  if (length(fasta.files) == 0){
    system(paste0("rm -rf ", irma.directory, "/", sample.names[i]))
    next
  } 
  
  # Check whether the file exists and overwrite it when requested.
  if (file.exists(paste0(output.directory, "/", sample.names[i], ".fasta")) == TRUE){
    system(paste0("rm ", output.directory, "/", sample.names[i], ".fasta"))
  }
  
  # Concatenate all FASTA files.
  system(paste0(
    "cat ", paste0(fasta.files, collapse = " "),
    " > ", output.directory, "/", sample.names[i], ".fasta"
  ))

  #################################################
  ### Part B: collect the amended consensus
  #################################################
  # IRMA writes the amended consensus as *.fa in amended_consensus/: padded,
  # low-coverage positions masked to N, minor alleles as IUPAC codes. Collect it
  # into IRMA-amended-contigs/.
  amended.files = list.files(
    paste0(irma.directory, "/", sample.names[i], "/amended_consensus"),
    pattern = "\\.(fa|fasta)$", full.names = T
  )

  if (length(amended.files) != 0){
    if (file.exists(paste0(amended.directory, "/", sample.names[i], ".fasta")) == TRUE){
      system(paste0("rm ", amended.directory, "/", sample.names[i], ".fasta"))
    }
    system(paste0(
      "cat ", paste0(amended.files, collapse = " "),
      " > ", amended.directory, "/", sample.names[i], ".fasta"
    ))
  } # end amended block

}# end i loop



# End of script.
