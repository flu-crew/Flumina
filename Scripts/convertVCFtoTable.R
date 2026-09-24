#### Run the Flumina pipeline first to obtain variant calls in VCF files.

# Required packages:
# 1. data.table package

# Runtime is approximately 2–3 minutes for 500 samples.

args = commandArgs(trailingOnly = TRUE)
# args = "config.cfg"  # Use a local configuration file during development.

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

# Output directory for the analysis.
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/variant_analysis")

# VCF directory name; provide the full path when it is not in the working
# directory.
vcf.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/vcf_files")

# Output table name.
save.name = "variant-table"

# Amino-acid positions come from the actual coding intervals, not from
# ceiling(position/3) - see Scripts/fluORFs.R for why that is wrong for M2,
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

# Return the config value for a key, or the default when the key is missing,
# blank or NULL. findAAChanges.R uses the same helper.
cfg = function(key, default = NULL) {
  v = config[[key]]
  if (is.null(v)) return(default)
  v = gsub("\"", "", trimws(v))
  if (!nzchar(v) || v == "NULL") return(default)
  v
}

# Thresholds for the call assessment at the end of this script. They use the
# same config keys as FluLens loadRunThresholds. FluPore writes MIN_FREQ and
# Flumina writes MIN_ALLELE_FREQUENCY.
freq.key  = if (is.null(config[["MIN_FREQ"]])) "MIN_ALLELE_FREQUENCY" else "MIN_FREQ"
min.depth = suppressWarnings(as.numeric(cfg("MIN_DEPTH", "100")))
min.alt   = suppressWarnings(as.numeric(cfg("MIN_ALT", "10")))
min.freq  = suppressWarnings(as.numeric(cfg(freq.key, "0.01")))

# Stop when a threshold in the config is not a number, because the script cannot
# apply it.
if (is.na(min.depth)) { stop("MIN_DEPTH in the config is not a number: ", cfg("MIN_DEPTH")) }
if (is.na(min.alt)) { stop("MIN_ALT in the config is not a number: ", cfg("MIN_ALT")) }
if (is.na(min.freq)) { stop(freq.key, " in the config is not a number: ", cfg(freq.key)) }

#output.directory = "/Volumes/Extreme_SSD/Bailey_project/variant_analysis"
#vcf.directory = "/Volumes/Extreme_SSD/Bailey_project/vcf_files"


#############################################
#### Should not need to modify below here
#############################################

#### LoFreq

# VCF file name or path for analysis.
vcf.string = "lofreq-called-variants.vcf" #or "gatk4-filtered-snps.vcf"

# Create the output directory.
dir.create(output.directory)

# Collect the multi-file inputs.
all.files = list.files(vcf.directory, recursive = T)
vcf.files = all.files[grep(paste0(vcf.string, "$"), all.files)]

# Initialize the collected variant data.
header.data = c("method", "sample", "locus", "position", "reference",
                "alternative", "quality", "depth", "map_quality", "allele_frequency", "aa_position")

# Initialize the data frame used to collect records.
collect.data = data.table::data.table(matrix(as.numeric(0),
                                             nrow = length(vcf.files)*1000,
                                             ncol = length(header.data)))
data.table::setnames(collect.data, header.data)

collect.data[, method:=as.character(method)]
collect.data[, sample:=as.character(sample)]
collect.data[, locus:=as.character(locus)]
collect.data[, reference:=as.character(reference)]
collect.data[, alternative:=as.character(alternative)]

# Iterate over loci and process their records.
# seq_along, NOT 1:length(). With no matching VCFs length() is 0 and 1:0 is
# c(1, 0), so the loop RUNS, vcf.files[1] is NA, and the script dies with
# "cannot open file '<outdir>/vcf_files/NA': No such file or directory" - an
# error naming a file that was never meant to exist, which says nothing about
# the real problem, that the directory matched nothing.
x = 1
for (i in seq_along(vcf.files)){
  
  # Count comment lines to locate the first data line.
  VCF = file(paste0(vcf.directory, "/", vcf.files[i]), "r")
  skip = 0
  line = readLines(VCF, 1)
  
  # Locate the contig definition line.
  while(!grepl("#CHROM", line)) {
    skip = skip + 1
    line = readLines(VCF, 1)
  }
  
  close(VCF)
  
  # Read the VCF after determining how many lines to skip.
  VCF = read.table(paste0(vcf.directory, "/", vcf.files[i]), skip = skip, comment.char = "", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)
  
  if (is.null(nrow(VCF)) == TRUE || nrow(VCF) == 0){
    
    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "LoFreq")
    # Collect the variant data.
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    # Extract the sample data.
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = 0 )
    # Extract the sequence length.
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
    
    # Collect the variant data.
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    # Extract the sample data.
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = VCF$'#CHROM'[j] )
    # Extract the sequence length.
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = VCF$POS[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = VCF$REF[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = VCF$ALT[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = VCF$QUAL[j] )
    
    # Extract the read depth.
    depth.val = as.numeric(gsub(";", "", gsub(";.*", "", gsub(".*DP=", "", VCF[j,]$INFO))) ) 
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = depth.val)
    
    # Extract the mapping quality.
    if (length(grep("MQ=", VCF[j,]$INFO)) != 0){
      mq.val = as.numeric(gsub(";.*", "", gsub(".*;MQ=", "", VCF[j,]$INFO)))
    } else { mq.val = NA }
    
    # Extract the mapping quality.
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = mq.val)

    # Extract the allele frequency.
    freq.val = as.numeric(gsub(";.*", "", gsub(".*;AF=", "", VCF[j,]$INFO)))
    data.table::set(collect.data, i = as.integer(x), j = match("allele_frequency", header.data), value = freq.val)
    
    # aa_position is NOT computed here any more. It depends on which product a
    # position codes for, and a position can code for two, so it is filled in
    # after both callers have been read - see the annotation block below.
    # Increment the record counter.
    x = x + 1
  } # End of j loop.
  
} # End of i loop.

# Remove samples with no records.
collect.data = collect.data[collect.data$sample != 0,]
collect.data = collect.data[collect.data$locus != 0,]

# Prevent amino acid T from being converted to TRUE.
collect.data$alternative[collect.data$alternative == "TRUE"] = "T"
collect.data$reference[collect.data$reference == "TRUE"] = "T"

lofreq.data = collect.data

#############################################
#### iVar
#############################################
# iVar writes a TAB-SEPARATED TABLE, not a VCF, so none of the VCF machinery
# above applies - there is no #CHROM line to skip to and no INFO field to
# regex apart. It is read directly.
#
# The file may legitimately be absent: IVAR=FALSE means the process never ran,
# and this has to degrade to "no iVar rows" rather than fail, the same way an
# absent metadata file does.
#
# Everything is read as CHARACTER and converted explicitly. read.table's type
# guessing is what turns a REF/ALT column of A/C/G/T into logicals when the
# only values present are T - the same "amino acid T becomes TRUE" problem the
# VCF blocks patch up afterwards. Reading as character avoids it at the source,
# and iVar's own PASS column is genuinely TRUE/FALSE, so a column of "T" and a
# column of booleans would otherwise be indistinguishable to the reader.
#
# iVar's columns beyond PASS (GFF_FEATURE, REF_CODON, REF_AA, ALT_CODON,
# ALT_AA, POS_AA) are always present and always NA here, because the process
# deliberately passes no GFF: iVar translates from the start of each reference
# sequence in frame 1, which is wrong for M2, NEP, PA-X and PB1-F2 in exactly
# the way ceiling(POS/3) was. The amino-acid annotation below is authoritative
# for every caller.

ivar.string = "ivar-called-variants.tsv"
all.ivar.files = list.files(vcf.directory, recursive = TRUE)
ivar.files = all.ivar.files[grep(paste0(ivar.string, "$"), all.ivar.files)]

ivar.rows   = list()
n.ivar.indel = 0
# Indel positions for the dist_to_indel column, and the samples an iVar file was
# actually found for. The second is what separates "no indel near this call"
# from "no way to tell" - see the indel-proximity block after the callers merge.
indel.pos = list()
ivar.samples.seen = character(0)

for (i in seq_along(ivar.files)) {

  ivar.path = paste0(vcf.directory, "/", ivar.files[i])

  # iVar TSVs are occasionally column-RAGGED: a few indel/edge rows carry an
  # extra tab-delimited field (e.g. 21 columns against the 20-column header, with
  # an empty REF/ALT). A plain read.table then throws "more columns than column
  # names"; the try() below swallows it and the ENTIRE sample's calls are dropped,
  # turning a real run into an empty variant table. Read the lines first and keep
  # only those whose field count matches the header - the discarded rows are
  # malformed indel artefacts the indel filter would drop anyway. Count by tab
  # (robust to trailing empty fields, which strsplit would silently eat).
  raw = readLines(ivar.path, warn = FALSE)
  if (length(raw) < 1L) { next }
  n.fields = function(s) vapply(gregexpr("\t", s, fixed = TRUE),
                                function(g) if (g[1L] == -1L) 1L else length(g) + 1L,
                                integer(1L))
  ncol.header = n.fields(raw[1L])
  ivar.body   = raw[-1L]
  keep.rows   = n.fields(ivar.body) == ncol.header
  if (any(!keep.rows)) {
    cat(sprintf("  note: %s -- dropped %d ragged row(s) (fields != %d)\n",
                basename(ivar.path), sum(!keep.rows), ncol.header))
  }
  ivar.tab = try(utils::read.table(text = c(raw[1L], ivar.body[keep.rows]),
                                   sep = "\t", header = TRUE,
                                   stringsAsFactors = FALSE, check.names = FALSE,
                                   colClasses = "character", quote = "",
                                   comment.char = "", na.strings = ""),
                 silent = TRUE)

  if (inherits(ivar.tab, "try-error") || is.null(nrow(ivar.tab)) || nrow(ivar.tab) == 0) { next }

  ivar.sample = gsub("/.*", "", ivar.files[i])

  # iVar reports indels in ALT as "+A" / "-T". Every downstream coordinate here
  # is SNP-based - an indel would shift the reading frame of everything after it
  # and there is no alignment to shift against - so they are dropped and
  # counted, never silently bounds-checked away. The other two callers' indels
  # are separated upstream (SelectVariants / filter_INDEL); iVar puts both kinds
  # in one file, so this is where it happens for iVar.
  is.indel = grepl("^[+-]", ivar.tab$ALT)
  n.ivar.indel = n.ivar.indel + sum(is.indel)

  # Their POSITIONS are kept even though the rows are dropped: `lofreq call`
  # runs with -B, and the calls that buys are enriched for sitting next to
  # indels, so every call carries its distance to the nearest one.
  #
  # iVar is the ONLY usable source - GATK4 is a genotype caller and misses
  # indels below genotype frequency, so its indel VCF is near-empty and would
  # report that nothing anywhere is indel-adjacent. Numbers in HANDOFF.md.
  #
  # Recorded before ivar.tab is subset and before the `next` below, because a
  # sample whose iVar output is ALL indels still has indels.
  ivar.samples.seen = c(ivar.samples.seen, ivar.sample)
  if (any(is.indel)) {
    indel.pos[[length(indel.pos) + 1]] = data.frame(
      sample   = ivar.sample,
      locus    = as.character(ivar.tab$REGION[is.indel]),
      position = as.numeric(ivar.tab$POS[is.indel]),
      stringsAsFactors = FALSE)
  }

  ivar.tab = ivar.tab[!is.indel, , drop = FALSE]

  if (nrow(ivar.tab) == 0) { next }

  ivar.rows[[length(ivar.rows) + 1]] = data.frame(
    method           = "iVar",
    sample           = ivar.sample,
    locus            = ivar.tab$REGION,
    position         = as.numeric(ivar.tab$POS),
    reference        = ivar.tab$REF,
    alternative      = ivar.tab$ALT,
    # ALT_QUAL is the mean quality of the reads carrying the alt, which is the
    # closest thing iVar reports to the VCF QUAL the other two supply.
    quality          = as.numeric(ivar.tab$ALT_QUAL),
    depth            = as.numeric(ivar.tab$TOTAL_DP),
    # iVar reports no mapping quality. NA, not 0 - 0 would read as "measured and
    # terrible" rather than "not measured", and MIN_QUALITY filters on it.
    map_quality      = NA_real_,
    allele_frequency = as.numeric(ivar.tab$ALT_FREQ),
    aa_position      = 0,
    # iVar's own verdict on its own call, from its Fisher exact test against the
    # sequencing error rate. Carried because it is the same situation as GATK4's
    # FILTER column: iVar ANNOTATES rather than removes, so without this a flagged
    # row arrives looking like a clean call.
    ivar_pass        = ivar.tab$PASS,
    # GATK4's FILTER verdict, which iVar rows do not have. Declared here rather
    # than assigned after the bind because ivar.data is legitimately empty when
    # IVAR=FALSE, and adding a column to a zero-row table is the kind of thing
    # that works until the day it is actually zero-row.
    gatk_filter      = NA_character_,
    stringsAsFactors = FALSE
  )
}

if (length(ivar.rows) > 0) {
  ivar.data = data.table::rbindlist(ivar.rows)
} else {
  # Same columns, no rows, so the rbind below behaves identically whether iVar
  # ran or not.
  ivar.data = data.table::data.table(
    method = character(), sample = character(), locus = character(),
    position = numeric(), reference = character(), alternative = character(),
    quality = numeric(), depth = numeric(), map_quality = numeric(),
    allele_frequency = numeric(), aa_position = numeric(), ivar_pass = character(),
    gatk_filter = character())
}

cat(sprintf("iVar: %d files, %d SNP rows kept, %d indel rows dropped\n",
            length(ivar.files), nrow(ivar.data), n.ivar.indel))

#############################################
#### GATK4
#############################################

# VCF file name or path for analysis.
vcf.string = "gatk4-filtered-snps.vcf" #or "gatk4-filtered-snps.vcf"

# Create the output directory.
dir.create(output.directory)

# Collect the multi-file inputs.
all.files = list.files(vcf.directory, recursive = T)
vcf.files = all.files[grep(paste0(vcf.string, "$"), all.files)]

# Initialize the collected variant data.
# gatk_filter carries VariantFiltration's own verdict - see the block where it
# is read, below.
header.data = c("method", "sample", "locus", "position", "reference",
                "alternative", "quality", "depth", "map_quality", "allele_frequency", "aa_position",
                "gatk_filter")

# Initialize the data frame used to collect records.
collect.data = data.table::data.table(matrix(as.numeric(0),
                                             nrow = length(vcf.files)*1000,
                                             ncol = length(header.data)))
data.table::setnames(collect.data, header.data)

collect.data[, method:=as.character(method)]
collect.data[, sample:=as.character(sample)]
collect.data[, locus:=as.character(locus)]
collect.data[, reference:=as.character(reference)]
collect.data[, alternative:=as.character(alternative)]
collect.data[, gatk_filter:=as.character(gatk_filter)]

# Iterate over loci and process their records.
# seq_along, NOT 1:length(). With no matching VCFs length() is 0 and 1:0 is
# c(1, 0), so the loop RUNS, vcf.files[1] is NA, and the script dies with
# "cannot open file '<outdir>/vcf_files/NA': No such file or directory" - an
# error naming a file that was never meant to exist, which says nothing about
# the real problem, that the directory matched nothing.
x = 1
for (i in seq_along(vcf.files)){

  # Count comment lines to locate the first data line.
  VCF = file(paste0(vcf.directory, "/", vcf.files[i]), "r")
  skip = 0
  line = readLines(VCF, 1)

  # Locate the contig definition line.
  while(!grepl("#CHROM", line)) {
    skip = skip + 1
    line = readLines(VCF, 1)
  }

  close(VCF)

  # Read the VCF after determining how many lines to skip.
  VCF = read.table(paste0(vcf.directory, "/", vcf.files[i]), skip = skip, comment.char = "", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)

  # FILTER is a mandatory VCF column, but read it defensively and once per file
  # rather than per row. Values are kept VERBATIM, including a bare "." - that
  # is the VCF's own way of saying "no filtering was applied", which is a
  # different statement from the NA the other two callers get, and collapsing
  # the two would lose exactly the distinction this column exists to make.
  filter.col = if ("FILTER" %in% names(VCF)) as.character(VCF$FILTER) else rep(NA_character_, max(nrow(VCF), 0))

  if (is.null(nrow(VCF)) == TRUE || nrow(VCF) == 0){

    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "GATK4")
    # Collect the variant data.
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    # Extract the sample data.
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = 0 )
    # Extract the sequence length.
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = 0 )
    data.table::set(collect.data, i = as.integer(x), j = match("aa_position", header.data), value = 0)
    data.table::set(collect.data, i = as.integer(x), j = match("gatk_filter", header.data), value = 0)
    x = x + 1
    next
  }
  
  for (j in 1:nrow(VCF)){
    
    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "GATK4")

    # GATK4's own verdict on its own call, and the reason it is worth a column
    # is the same one that earned ivar_pass its own: gatk4-filtered-snps.vcf is
    # VariantFiltration's output, which ANNOTATES rather than removes. Downstream
    # is expected to honour FILTER, so without this column a flagged record looks
    # exactly like a clean call.
    #
    # Annotated, not dropped. Whether non-PASS rows should be removed outright
    # is a decision that would change published results, and it is not this
    # script's to make - but it cannot be made at all while the column is
    # invisible.
    data.table::set(collect.data, i = as.integer(x), j = match("gatk_filter", header.data), value = filter.col[j])

    # Collect the variant data.
    data.table::set(collect.data, i = as.integer(x), j = match("sample", header.data), value = gsub("/.*", "", vcf.files[i]) )
    # Extract the sample data.
    data.table::set(collect.data, i = as.integer(x), j = match("locus", header.data), value = VCF$'#CHROM'[j] )
    # Extract the sequence length.
    data.table::set(collect.data, i = as.integer(x), j = match("position", header.data), value = VCF$POS[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("reference", header.data), value = VCF$REF[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("alternative", header.data), value = VCF$ALT[j] )
    data.table::set(collect.data, i = as.integer(x), j = match("quality", header.data), value = VCF$QUAL[j] )
    
    # Extract the read depth.
    depth.val = as.numeric(gsub(";", "", gsub(";.*", "", gsub(".*DP=", "", VCF[j,]$INFO))) ) 
    data.table::set(collect.data, i = as.integer(x), j = match("depth", header.data), value = depth.val)
    
    # Extract the mapping quality.
    if (length(grep("MQ=", VCF[j,]$INFO)) != 0){
      mq.val = as.numeric(gsub(";.*", "", gsub(".*;MQ=", "", VCF[j,]$INFO)))
    } else { mq.val = NA }
    
    # Extract the mapping quality.
    data.table::set(collect.data, i = as.integer(x), j = match("map_quality", header.data), value = mq.val)
    
    # Extract the allele frequency.
    freq.val = as.numeric(gsub(";.*", "", gsub(".*;AF=", "", VCF[j,]$INFO)))
    data.table::set(collect.data, i = as.integer(x), j = match("allele_frequency", header.data), value = freq.val)
    
    # aa_position is NOT computed here any more. It depends on which product a
    # position codes for, and a position can code for two, so it is filled in
    # after both callers have been read - see the annotation block below.
    # Increment the record counter.
    x = x + 1
  } # End of j loop.
  
} # End of i loop.

# Remove samples with no records.
collect.data = collect.data[collect.data$sample != 0,]
collect.data = collect.data[collect.data$locus != 0,]

# Prevent amino acid T from being converted to TRUE.
collect.data$alternative[collect.data$alternative == "TRUE"] = "T"
collect.data$reference[collect.data$reference == "TRUE"] = "T"

# ivar_pass belongs only to iVar rows and gatk_filter only to GATK4's, so each
# caller is given the other's column as NA before the bind. NA means "this
# caller has no such verdict" - a different statement from iVar's FALSE, and a
# different one again from GATK4's own ".", which means "no filter was
# applied". None of the three may collapse into another.
#
# set() rather than $<-, because it adds the column by reference and is correct
# on a zero-row table; $<- on a zero-row data.table is an error waiting for the
# first run that produces no calls for a caller.
data.table::set(lofreq.data,  j = "ivar_pass",   value = NA_character_)
data.table::set(lofreq.data,  j = "gatk_filter", value = NA_character_)
data.table::set(collect.data, j = "ivar_pass",   value = NA_character_)

# rbind binds by POSITION, not by name. gatk_filter arrives inside GATK4's own
# header.data but is appended to the end of the other two, so without this the
# three tables agree on their column NAMES and disagree on their order - which
# silently interleaves ivar_pass and gatk_filter values between callers.
bind.order = c("method", "sample", "locus", "position", "reference", "alternative",
               "quality", "depth", "map_quality", "allele_frequency", "aa_position",
               "ivar_pass", "gatk_filter")
data.table::setcolorder(lofreq.data,  bind.order)
data.table::setcolorder(ivar.data,    bind.order)
data.table::setcolorder(collect.data, bind.order)

final.data = rbind(lofreq.data, ivar.data, collect.data)

#############################################
#### Reconciling the two callers
#############################################
# LoFreq's allele_frequency is an allele FRACTION. GATK4's is a GENOTYPE - a
# hom-alt call is 1.0 whatever the reads say. LoFreq tracks the observed read
# fraction; GATK4 does not report the same quantity.
#
# Both callers also emit a row for the same change, so at most sites the table
# carries TWO ROWS PER VARIANT. Counting rows therefore over-counts changes: a
# single change can appear twice.
#
# Nothing existing is rewritten. `allele_frequency` keeps exactly what each
# caller reported, because silently changing what a published column means is
# worse than the situation it would fix. Three columns are ADDED so that every
# downstream consumer stops having to work this out for itself:
#
#   variant_id       identical for every row describing the same change
#   af_type          "fraction" or "genotype" - what allele_frequency IS
#   allele_fraction  the best available true fraction: LoFreq's own value, or
#                    LoFreq's borrowed for a matching GATK4 row, NA when only
#                    GATK4 saw it and there is nothing to borrow
#
# Run BEFORE the amino-acid annotation below, which expands to one row per
# product. After that expansion a variant_id legitimately appears once per
# product (the same nucleotide read in two frames), so dedupe on
# variant_id + product, never variant_id alone.

final.data$variant_id = paste(final.data$sample, final.data$locus,
                              final.data$position, final.data$alternative,
                              sep = "|")
# iVar reports an allele FRACTION, like LoFreq and unlike GATK4. So there are
# now two callers measuring the same quantity and one measuring a different one,
# and af_type says which - never the caller name, which is what a consumer would
# otherwise have to hardcode a list against.
final.data$af_type = ifelse(final.data$method == "GATK4", "genotype", "fraction")

is.lofreq = final.data$method == "LoFreq"
is.ivar   = final.data$method == "iVar"

lofreq.af = stats::setNames(as.numeric(final.data$allele_frequency[is.lofreq]),
                            final.data$variant_id[is.lofreq])
ivar.af   = stats::setNames(as.numeric(final.data$allele_frequency[is.ivar]),
                            final.data$variant_id[is.ivar])

# LoFreq keeps priority for allele_fraction. Its value matches the observed read
# fraction, and it is what existing analyses of these tables were built on. iVar
# fills the slot only where LoFreq never saw the change.
borrowed.af = unname(lofreq.af[final.data$variant_id])
borrowed.af = ifelse(is.na(borrowed.af), unname(ivar.af[final.data$variant_id]), borrowed.af)

final.data$allele_fraction = ifelse(is.lofreq | is.ivar,
                                    as.numeric(final.data$allele_frequency),
                                    borrowed.af)

# Where BOTH fraction callers saw the same change, record how far apart they
# are. Not a boolean and not a correction: the two are independent measurements
# of one quantity, so the gap between them is evidence about the call, and
# collapsing it to a flag at a threshold chosen here would decide for every
# downstream consumer what "disagreement" means. NA where only one caller saw
# the change, which is different from agreeing.
lo.for.row = unname(lofreq.af[final.data$variant_id])
iv.for.row = unname(ivar.af[final.data$variant_id])
final.data$af_conflict = ifelse(is.na(lo.for.row) | is.na(iv.for.row),
                                NA_real_, abs(lo.for.row - iv.for.row))

n.changes = length(unique(final.data$variant_id))
n.dup     = sum(duplicated(final.data$variant_id))
n.borrow  = sum(!is.lofreq & !is.ivar & !is.na(final.data$allele_fraction))
n.geno    = sum(!is.lofreq & !is.ivar &  is.na(final.data$allele_fraction))
cat(sprintf("Caller reconciliation: %d rows -> %d distinct changes (%d duplicate reports)\n",
            nrow(final.data), n.changes, n.dup))
cat(sprintf("  Rows by caller: LoFreq %d, iVar %d, GATK4 %d\n",
            sum(is.lofreq), sum(is.ivar), sum(!is.lofreq & !is.ivar)))
cat(sprintf("  GATK4 rows: %d took a fraction from LoFreq or iVar, %d genotype-only (allele_fraction NA)\n",
            n.borrow, n.geno))

# Counted per CALL, before the amino-acid expansion below, so this is a count of
# records GATK4 wrote rather than of table rows. Print the FILTER breakdown so a
# heavily flagged run is visible.
gatk.filt = final.data$gatk_filter[final.data$method == "GATK4"]
if (length(gatk.filt) > 0) {
  ft = sort(table(ifelse(is.na(gatk.filt), "(missing)", gatk.filt)), decreasing = TRUE)
  cat(sprintf("  GATK4 FILTER: %s\n", paste(names(ft), ft, sep = "=", collapse = "  ")))
  n.pass = sum(gatk.filt %in% c("PASS", "."), na.rm = TRUE)
  cat(sprintf("    %d of %d GATK4 records PASS or unfiltered, %d flagged by VariantFiltration\n",
              n.pass, length(gatk.filt), length(gatk.filt) - n.pass))
}

# Reported per CHANGE rather than per row: both callers contribute a row each,
# so counting rows would double every disagreement.
conflict.by.id = final.data$af_conflict[!duplicated(final.data$variant_id)]
n.shared = sum(!is.na(conflict.by.id))
if (n.shared > 0) {
  cat(sprintf("  LoFreq vs iVar: %d changes seen by both, median |difference| %.4f, max %.4f\n",
              n.shared, stats::median(conflict.by.id, na.rm = TRUE),
              max(conflict.by.id, na.rm = TRUE)))
  cat(sprintf("    %d differ by >1 percentage point, %d by >5\n",
              sum(conflict.by.id > 0.01, na.rm = TRUE),
              sum(conflict.by.id > 0.05, na.rm = TRUE)))
}

#############################################
#### Amino-acid annotation, through the real coding intervals
#############################################
# Eight segments, twelve proteins. A nucleotide inside PA may also code for
# PA-X in a different frame; one inside MP past nucleotide 756 codes for M2 and
# for nothing else. So this expands to ONE ROW PER PRODUCT rather than assuming
# a single reading frame starting at nucleotide 1.
#
# `locus` deliberately stays the SEGMENT. Everything downstream that keys on it
# - the curated-database join in outputSummary.R above all - keeps working
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

#############################################
#### Indel proximity
#############################################
# `lofreq call` runs with -B (BAQ off), which recovers real calls but also
# admits ones next to indels - what BAQ was suppressing. The cost is carried as
# a column rather than paid in lost calls, the same way gatk_filter and
# ivar_pass annotate rather than remove.
#
# TWO columns, because one cannot say both things:
#   dist_to_indel  bases to the nearest iVar indel in the same sample+locus,
#                  NA when that segment has none
#   indel_source   whether an iVar file existed for that sample at all
#
# Without the second, NA is ambiguous between "measured, nothing near it" and
# "iVar never ran, so nobody looked" - opposite conclusions. Same argument as
# min_depth riding on every row of the IRMA frame table.
#
# Computed AFTER the amino-acid annotation: distance is a property of the
# nucleotide position, so it survives that step's one-row-per-ORF expansion,
# and this way it does not rely on flu_annotate_positions carrying columns it
# knows nothing about.
indel.table = if (length(indel.pos) > 0) do.call(rbind, indel.pos) else NULL

dist.vec = rep(NA_real_, nrow(final.data))
if (!is.null(indel.table) && nrow(indel.table) > 0) {
  # "\r" as the key separator: sample and locus names cannot contain it, so no
  # pair of real names can collide into one key.
  indel.key = paste(indel.table$sample, indel.table$locus, sep = "\r")
  by.key    = split(indel.table$position, indel.key)
  row.key   = paste(final.data$sample, final.data$locus, sep = "\r")
  hit       = which(row.key %in% names(by.key))
  if (length(hit) > 0) {
    dist.vec[hit] = vapply(hit, function(j)
      min(abs(final.data$position[j] - by.key[[row.key[j]]])), numeric(1))
  }
}

final.data$dist_to_indel = dist.vec
final.data$indel_source  = final.data$sample %in% ivar.samples.seen

n.src   = sum(final.data$indel_source)
n.near  = sum(!is.na(dist.vec) & dist.vec <= 10)
cat(sprintf("Indel proximity: %d indel position(s) from %d sample(s) with an iVar source; %d of %d rows within 10 bp of one\n",
            if (is.null(indel.table)) 0L else nrow(indel.table),
            length(unique(ivar.samples.seen)), n.near, nrow(final.data)))
if (n.src < nrow(final.data))
  cat(sprintf("  %d row(s) have NO indel source (iVar did not run for that sample): dist_to_indel is NA and means UNKNOWN, not far\n",
              nrow(final.data) - n.src))
#############################################
#### Call assessment
#############################################
# Give each call the same verdict as FluLens assessCore. The verdict is set here,
# not in FluLens, so every tool that reads the table gets the same verdict.
# The verdict uses the DP4 read counts of the call: reference forward, reference
# reverse, alt forward and alt reverse. DP4 is not a column, because its commas
# break this unquoted CSV. min.depth, min.alt and min.freq come from the config at
# the top of this script. This block adds three columns:
#   alt_reads     reads that support the alt (DP4 alt forward + alt reverse)
#   strand_class  balanced, some-skew, skewed, too-few-alt, no-ref-control,
#                 not-assessed (strand test off) or NA (no DP4 record)
#   assessment    Looks real, Treat with caution, Likely artefact or Cannot assess

# Fixed strand limits. These are AS_STRAND_MIN_ALT, AS_SKEW_BAD and AS_SKEW_WARN
# in FluLens.
strand.min.reads = 4
skew.bad         = 0.40
skew.warn        = 0.25

# ONT reads have a strand bias that is not an artefact, so FluPore runs do not
# use the strand test. A MIN_FREQ key identifies a FluPore run, as in FluLens.
run.strand.bias = is.null(config[["MIN_FREQ"]])

# Read the DP4 counts from the same per-sample files that FluLens reads. The id
# of a record is the caller, sample, locus and position. The iVar id also has
# the alt base, because iVar writes one row for each alt base.
dp4.rows = list()

# LoFreq writes DP4 in the INFO field of each VCF record.
dp4.pattern = ".*DP4=([0-9]+,[0-9]+,[0-9]+,[0-9]+).*"
lofreq.files = list.files(vcf.directory, pattern = "lofreq-called-variants.vcf$", recursive = TRUE)
for (i in seq_along(lofreq.files)) {
  vcf.lines = readLines(paste0(vcf.directory, "/", lofreq.files[i]), warn = FALSE)
  vcf.lines = vcf.lines[!startsWith(vcf.lines, "#") & grepl(dp4.pattern, vcf.lines)]
  if (length(vcf.lines) == 0) { next }

  # Split the records into VCF columns: 1 is CHROM, 2 is POS and 8 is INFO.
  fields = data.table::tstrsplit(vcf.lines, "\t", fixed = TRUE)
  dp4 = sub(dp4.pattern, "\\1", fields[[8]])
  dp4.counts = data.table::tstrsplit(dp4, ",", fixed = TRUE)

  dp4.rows[[length(dp4.rows) + 1]] = data.frame(
    id      = paste("LoFreq", gsub("/.*", "", lofreq.files[i]), fields[[1]],
                    as.numeric(fields[[2]]), sep = "|"),
    ref.fwd = as.numeric(dp4.counts[[1]]),
    ref.rev = as.numeric(dp4.counts[[2]]),
    alt.fwd = as.numeric(dp4.counts[[3]]),
    alt.rev = as.numeric(dp4.counts[[4]]),
    stringsAsFactors = FALSE)
}

# iVar writes the counts in the columns of its TSV. ivar.files is the file list
# from the iVar block above.
ivar.need = c("REGION", "POS", "ALT", "REF_DP", "REF_RV", "ALT_DP", "ALT_RV")
for (i in seq_along(ivar.files)) {
  ivar.tab = try(data.table::fread(paste0(vcf.directory, "/", ivar.files[i]), sep = "\t", fill = TRUE,
                                   colClasses = list(character = c("REGION", "REF", "ALT"))),
                 silent = TRUE)
  if (inherits(ivar.tab, "try-error") || nrow(ivar.tab) == 0 ||
      !all(ivar.need %in% names(ivar.tab))) { next }

  # The table has no indels, so skip them here too.
  ivar.tab = ivar.tab[!grepl("^[+-]", ivar.tab$ALT), ]
  if (nrow(ivar.tab) == 0) { next }

  # iVar gives the total and the reverse count, so forward = total - reverse.
  dp4.rows[[length(dp4.rows) + 1]] = data.frame(
    id      = paste("iVar", gsub("/.*", "", ivar.files[i]), ivar.tab$REGION,
                    as.numeric(ivar.tab$POS), ivar.tab$ALT, sep = "|"),
    ref.fwd = as.numeric(ivar.tab$REF_DP) - as.numeric(ivar.tab$REF_RV),
    ref.rev = as.numeric(ivar.tab$REF_RV),
    alt.fwd = as.numeric(ivar.tab$ALT_DP) - as.numeric(ivar.tab$ALT_RV),
    alt.rev = as.numeric(ivar.tab$ALT_RV),
    stringsAsFactors = FALSE)
}

if (length(dp4.rows) > 0) {
  dp4.table = data.table::rbindlist(dp4.rows)
} else {
  # Same columns, no rows, so the lookup below works when no file has DP4.
  dp4.table = data.table::data.table(id = character(), ref.fwd = numeric(),
                                     ref.rev = numeric(), alt.fwd = numeric(),
                                     alt.rev = numeric())
}

# FluLens keeps the last record for an id, so do the same here.
dp4.table = dp4.table[!duplicated(dp4.table$id, fromLast = TRUE), ]

# Find the DP4 record of each call. iVar rows use the iVar record. LoFreq and
# GATK4 rows use the LoFreq record at the same position, because GATK4 writes no
# DP4 (as in FluLens strandRecOf). A call with no record gets NA counts.
row.id = ifelse(final.data$method == "iVar",
                paste("iVar", final.data$sample, final.data$locus,
                      final.data$position, final.data$alternative, sep = "|"),
                paste("LoFreq", final.data$sample, final.data$locus,
                      final.data$position, sep = "|"))
row.dp4 = match(row.id, dp4.table$id)
has.dp4 = !is.na(row.dp4)

ref.fwd   = dp4.table$ref.fwd[row.dp4]
alt.fwd   = dp4.table$alt.fwd[row.dp4]
ref.reads = ref.fwd + dp4.table$ref.rev[row.dp4]
alt.reads = alt.fwd + dp4.table$alt.rev[row.dp4]

# The forward share of the alt reads and of the reference reads. The skew is the
# difference between the two shares.
alt.frac = ifelse(alt.reads > 0, alt.fwd / alt.reads, 0)
ref.frac = ifelse(ref.reads > 0, ref.fwd / ref.reads, 0)
skew     = abs(alt.frac - ref.frac)

# Use allele_fraction. A GATK4 row with no fraction uses its allele_frequency.
freq = ifelse(is.na(final.data$allele_fraction), final.data$allele_frequency,
              final.data$allele_fraction)

# The strand class of each call. The tests are in the same order as in FluLens
# assessCore, and the first test that is true gives the class. A call with too
# few reference reads (a fixed call) gets no-ref-control, because the skew test
# needs reference reads.
strand.class = rep("not-assessed", nrow(final.data))
if (run.strand.bias) {
  strand.class = ifelse(alt.reads < strand.min.reads, "too-few-alt",
                 ifelse(alt.frac == 0 | alt.frac == 1, "skewed",
                 ifelse(ref.reads < strand.min.reads, "no-ref-control",
                 ifelse(skew > skew.bad, "skewed",
                 ifelse(skew > skew.warn, "some-skew", "balanced")))))
}
strand.class[!has.dp4] = NA_character_

# The verdict of each call, with the same tests as FluLens assessCore. A call
# with no DP4 record cannot be assessed.
likely.artefact = has.dp4 &
  (alt.reads < min.alt | (run.strand.bias & strand.class == "skewed"))
caution = has.dp4 & !likely.artefact &
  (final.data$depth < min.depth | freq < min.freq |
   (run.strand.bias & strand.class %in% c("some-skew", "too-few-alt")))
looks.real = has.dp4 & !likely.artefact & !caution

assessment = rep("Cannot assess", nrow(final.data))
assessment[likely.artefact] = "Likely artefact"
assessment[caution]         = "Treat with caution"
assessment[looks.real]      = "Looks real"

final.data$alt_reads    = alt.reads
final.data$strand_class = strand.class
final.data$assessment   = assessment

assessment.counts = table(factor(assessment,
                                 levels = c("Looks real", "Treat with caution",
                                            "Likely artefact", "Cannot assess")))
cat(sprintf("Call assessment (MIN_DEPTH=%g MIN_ALT=%g MIN_FREQ=%g strand_bias=%s): %s\n",
            min.depth, min.alt, min.freq, run.strand.bias,
            paste(names(assessment.counts), assessment.counts, sep = "=", collapse = "  ")))

# Save the data.
write.csv(final.data, paste0(output.directory, "/", save.name, ".csv"),
          row.names = F, quote = F)




#########################
###### End of script
#########################











