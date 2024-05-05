

#setwd("/Users/chutter/Dropbox/Research/0_Github/Flumina/test_example/")

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
read.directory = gsub("\"", "", config$READ_DIRECTORY)
rename.file = gsub("\"", "", config$RENAME_FILE)
overwrite = is.logical(gsub("\"", "", config$OVERWRITE))
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/organized-reads")

### Script start
#######################3

#Quick checks
options(stringsAsFactors = FALSE)
if (is.null(read.directory) == TRUE){ stop("Please provide a directory of raw reads.") }
if (file.exists(read.directory) == F){ stop("Input reads not found.") }
if (is.null(rename.file) == TRUE){ stop("Please provide a table of file to sample name conversions.") }

#Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
  if (overwrite == TRUE){
    system(paste0("rm -r ", output.directory))
    dir.create(output.directory)
  }
}#end else

#Read in sample data and finds reads
reads = list.files(read.directory, recursive = T, full.names = T)
reads = reads[grep("fastq.gz$|fastq$|fq.gz$|fq$", reads)]

sample.data = read.csv(rename.file)
if (nrow(sample.data) == 0){ return("no samples available to organize.") }

#Skips samples already finished
if (overwrite == FALSE){
  done.names = list.files(output.directory)
  sample.data = sample.data[!sample.data$Sample %in% done.names,]
} else { sample.data = sample.data }

if (nrow(sample.data) == 0){ return("no samples available to organize.") }

sample.names = unique(sample.data$Sample)

for (i in seq_along(sample.names)){

  temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]

  for (j in 1:nrow(temp.data)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    # Finds all files for this given sample and turns into a cat string
    sample.reads = reads[grep(pattern = paste0(temp.data$File[j], "_"), x = reads)]
    # Checks the Sample column in case already renamed
    if (length(sample.reads) == 0) {
      sample.reads = reads[grep(pattern = paste0(temp.data$Sample[j], "_"), x = reads)]
    }
    if (length(sample.reads) == 0) {
      sample.reads = reads[grep(pattern = temp.data$File[j], x = reads)]
    }
    if (length(sample.reads) == 0) {
      sample.reads = reads[grep(pattern = temp.data$Sample[j], x = reads)]
    }

    # Returns an error if reads are not found
    if (length(sample.reads) == 0) {
      print(paste0(
        temp.data$Sample[j], " does not have any reads present for files ",
        temp.data$File[j], " from the input spreadsheet."
      ))
	next
    } # end if statement

    if (length(sample.reads) == 1) {
      print(paste0(temp.data$Sample[j], " only one read file found. The other is missing."))
	next
    }

    if (length(sample.reads) >= 3) {
      print(paste0(temp.data$Sample[j], " had more than 2 read files associated, all file names must be unique."))
	next
    }

    #################################################
    ### Part B: Create directories and move files
    #################################################
    # Create sample directory
    out.path = paste0(output.directory, "/", temp.data$Sample[j])
    if (file.exists(out.path) == FALSE) {
      dir.create(out.path)
    }

    # sets up output reads
    outread.1 = paste0(out.path, "/", temp.data$Sample[j], "_L00", j, "_READ1.fastq.gz")
    outread.2 = paste0(out.path, "/", temp.data$Sample[j], "_L00", j, "_READ2.fastq.gz")

    system(paste0("cp ", sample.reads[1], " ", outread.1))
    system(paste0("cp ", sample.reads[2], " ", outread.2))
  } # end j loop

}#end i loop

