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
# seq_along, NOT 1:length(). With no matching VCFs length() is 0 and 1:0 is
# c(1, 0), so the loop RUNS, vcf.files[1] is NA, and the script dies with
# "cannot open file '<outdir>/vcf_files/NA': No such file or directory" — an
# error naming a file that was never meant to exist, which says nothing about
# the real problem, that the directory matched nothing.
x = 1
for (i in seq_along(vcf.files)){
  
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
#### iVar
#############################################
# iVar writes a TAB-SEPARATED TABLE, not a VCF, so none of the VCF machinery
# above applies — there is no #CHROM line to skip to and no INFO field to
# regex apart. It is read directly.
#
# The file may legitimately be absent: IVAR=FALSE means the process never ran,
# and this has to degrade to "no iVar rows" rather than fail, the same way an
# absent metadata file does.
#
# Everything is read as CHARACTER and converted explicitly. read.table's type
# guessing is what turns a REF/ALT column of A/C/G/T into logicals when the
# only values present are T — the same "amino acid T becomes TRUE" problem the
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
# from "no way to tell" — see the indel-proximity block after the callers merge.
indel.pos = list()
ivar.samples.seen = character(0)

for (i in seq_along(ivar.files)) {

  ivar.path = paste0(vcf.directory, "/", ivar.files[i])

  # iVar TSVs are occasionally column-RAGGED: a few indel/edge rows carry an
  # extra tab-delimited field (e.g. 21 columns against the 20-column header, with
  # an empty REF/ALT). A plain read.table then throws "more columns than column
  # names"; the try() below swallows it and the ENTIRE sample's calls are dropped,
  # turning a real run into an empty variant table. Read the lines first and keep
  # only those whose field count matches the header — the discarded rows are
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
  # is SNP-based — an indel would shift the reading frame of everything after it
  # and there is no alignment to shift against — so they are dropped and
  # counted, never silently bounds-checked away. The other two callers' indels
  # are separated upstream (SelectVariants / filter_INDEL); iVar puts both kinds
  # in one file, so this is where it happens for iVar.
  is.indel = grepl("^[+-]", ivar.tab$ALT)
  n.ivar.indel = n.ivar.indel + sum(is.indel)

  # Their POSITIONS are kept even though the rows are dropped: `lofreq call`
  # runs with -B, and the calls that buys are enriched for sitting next to
  # indels, so every call carries its distance to the nearest one.
  #
  # iVar is the ONLY usable source — GATK4 is a genotype caller and misses
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
    # iVar reports no mapping quality. NA, not 0 — 0 would read as "measured and
    # terrible" rather than "not measured", and MIN_QUALITY filters on it.
    map_quality      = NA_real_,
    allele_frequency = as.numeric(ivar.tab$ALT_FREQ),
    aa_position      = 0,
    # iVar's own verdict on its own call, from its Fisher exact test against the
    # sequencing error rate. Carried because it is the same situation as GATK4's
    # FILTER column: iVar ANNOTATES rather than removes, so without this every
    # row it flagged as failing arrives looking like a clean call. On MC-497's
    # real output, 25 of 34 rows are PASS=FALSE.
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

#the string or name of the VCF file for data anaylsis 
vcf.string = "gatk4-filtered-snps.vcf" #or "gatk4-filtered-snps.vcf"

#Creates output directory
dir.create(output.directory)

#Get multifile databases together
all.files = list.files(vcf.directory, recursive = T)
vcf.files = all.files[grep(paste0(vcf.string, "$"), all.files)]

#Collects the super cool data
# gatk_filter carries VariantFiltration's own verdict — see the block where it
# is read, below.
header.data = c("method", "sample", "locus", "position", "reference",
                "alternative", "quality", "depth", "map_quality", "allele_frequency", "aa_position",
                "gatk_filter")

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
collect.data[, gatk_filter:=as.character(gatk_filter)]

#Loops through each locus and does operations on them
# seq_along, NOT 1:length(). With no matching VCFs length() is 0 and 1:0 is
# c(1, 0), so the loop RUNS, vcf.files[1] is NA, and the script dies with
# "cannot open file '<outdir>/vcf_files/NA': No such file or directory" — an
# error naming a file that was never meant to exist, which says nothing about
# the real problem, that the directory matched nothing.
x = 1
for (i in seq_along(vcf.files)){

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

  # FILTER is a mandatory VCF column, but read it defensively and once per file
  # rather than per row. Values are kept VERBATIM, including a bare "." — that
  # is the VCF's own way of saying "no filtering was applied", which is a
  # different statement from the NA the other two callers get, and collapsing
  # the two would lose exactly the distinction this column exists to make.
  filter.col = if ("FILTER" %in% names(VCF)) as.character(VCF$FILTER) else rep(NA_character_, max(nrow(VCF), 0))

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
    data.table::set(collect.data, i = as.integer(x), j = match("gatk_filter", header.data), value = 0)
    x = x + 1
    next
  }
  
  for (j in 1:nrow(VCF)){
    
    data.table::set(collect.data, i = as.integer(x), j = match("method", header.data), value = "GATK4")

    # GATK4's own verdict on its own call, and the reason it is worth a column
    # is the same one that earned ivar_pass its own: gatk4-filtered-snps.vcf is
    # VariantFiltration's output, which ANNOTATES rather than removes.
    # Downstream is expected to honour FILTER, and until now nothing here read
    # it — so every record GATK4 flagged as failing its strand-bias tests
    # arrived in the table looking exactly like a clean call. On the swine WGS
    # run across all 143 samples that was EVERY record: 96 FS;SOR, 3 FS;QD;SOR,
    # 1 FS, and 0 PASS.
    #
    # Annotated, not dropped. Whether non-PASS rows should be removed outright
    # is a decision that would change published results, and it is not this
    # script's to make — but it cannot be made at all while the column is
    # invisible.
    data.table::set(collect.data, i = as.integer(x), j = match("gatk_filter", header.data), value = filter.col[j])

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

# ivar_pass belongs only to iVar rows and gatk_filter only to GATK4's, so each
# caller is given the other's column as NA before the bind. NA means "this
# caller has no such verdict" — a different statement from iVar's FALSE, and a
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
# three tables agree on their column NAMES and disagree on their order — which
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
# LoFreq's allele_frequency is an allele FRACTION. GATK4's is a GENOTYPE — a
# hom-alt call is 1.0 whatever the reads say. Measured against the BAM at
# PB2:1981 in MC-524 on the swine WGS run: LoFreq 85.98%, GATK4 100.00%, and
# 1,078 reads spanning the codon said 85.44%, with 12.89% still carrying the
# reference. LoFreq is within half a percent of the reads; GATK4 is not
# reporting the same quantity at all.
#
# Both callers also emit a row for the same change, so at most sites the table
# carries TWO ROWS PER VARIANT. Counting rows therefore over-counts changes: 88
# of the 381 codons on the swine run that looked like they carried multiple
# changes were one change reported twice.
#
# Nothing existing is rewritten. `allele_frequency` keeps exactly what each
# caller reported, because silently changing what a published column means is
# worse than the situation it would fix. Three columns are ADDED so that every
# downstream consumer stops having to work this out for itself:
#
#   variant_id       identical for every row describing the same change
#   af_type          "fraction" or "genotype" — what allele_frequency IS
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
# and af_type says which — never the caller name, which is what a consumer would
# otherwise have to hardcode a list against.
final.data$af_type = ifelse(final.data$method == "GATK4", "genotype", "fraction")

is.lofreq = final.data$method == "LoFreq"
is.ivar   = final.data$method == "iVar"

lofreq.af = stats::setNames(as.numeric(final.data$allele_frequency[is.lofreq]),
                            final.data$variant_id[is.lofreq])
ivar.af   = stats::setNames(as.numeric(final.data$allele_frequency[is.ivar]),
                            final.data$variant_id[is.ivar])

# LoFreq keeps priority for allele_fraction. It is the caller whose value was
# checked against the reads (85.98% against 85.44% measured at PB2:1981 in
# MC-524), it is what every existing analysis of these tables was built on, and
# a silent change of what this column means is worse than the situation it would
# fix. iVar fills the slot only where LoFreq never saw the change.
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
# records GATK4 wrote rather than of table rows. Printed because the swine WGS
# run's answer was that NOT ONE of 100 records passed — 96 FS;SOR, 3 FS;QD;SOR,
# 1 FS, 0 PASS — and that was invisible for as long as nothing read the column.
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

#############################################
#### Indel proximity
#############################################
# `lofreq call` runs with -B (BAQ off), which recovers real calls but also
# admits ones next to indels — what BAQ was suppressing. The cost is carried as
# a column rather than paid in lost calls, the same way gatk_filter and
# ivar_pass annotate rather than remove.
#
# TWO columns, because one cannot say both things:
#   dist_to_indel  bases to the nearest iVar indel in the same sample+locus,
#                  NA when that segment has none
#   indel_source   whether an iVar file existed for that sample at all
#
# Without the second, NA is ambiguous between "measured, nothing near it" and
# "iVar never ran, so nobody looked" — opposite conclusions. Same argument as
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
#### Call assessment — the FluLens verdict, stated in the table
#############################################
# The verdict is computed here, not in FluLens, so the table is the single source of
# truth and the viewer reads it. Same reason af_type / allele_fraction are stated
# above. The rule is FluLens assessCore, with one fix: the strand-skew test needs
# reference reads, so it is skipped when the reference side has fewer than
# AS_STRAND_MIN_ALT reads (a fixed call). FluLens uses the same rule.
# Three columns are added (raw DP4 is not, because its commas break this CSV):
#   alt_reads     reads that support the alt (DP4 alt fwd+rev), not depth*freq
#   strand_class  balanced / some-skew / skewed / too-few-alt / no-ref-control /
#                 not-assessed (ONT) / NA
#   assessment    Looks real / Treat with caution / Likely artefact / Cannot assess

# Thresholds, from the same config keys as FluLens loadRunThresholds.
num.cfg = function(key, default) {
  v = suppressWarnings(as.numeric(gsub("\"", "", config[[key]])))
  if (length(v) == 0L || is.na(v)) default else v
}
AS_MIN_DEPTH = num.cfg("MIN_DEPTH", 100)
AS_MIN_ALT   = num.cfg("MIN_ALT", 10)
# Flumina writes MIN_ALLELE_FREQUENCY; FluPore (ONT) writes MIN_FREQ. Match both.
AS_MIN_FREQ  = if (!is.null(config[["MIN_FREQ"]])) num.cfg("MIN_FREQ", 0.01) else num.cfg("MIN_ALLELE_FREQUENCY", 0.01)
AS_STRAND_MIN_ALT = 4L; AS_SKEW_BAD = 0.40; AS_SKEW_WARN = 0.25
# ONT reads have inherent strand bias, so the strand test is off for FluPore runs.
# MIN_FREQ marks one, as in FluLens runStrandBias.
run.strand.bias = is.null(config[["MIN_FREQ"]])

# DP4 for every call, from the same per-sample files FluLens reads.
read.lofreq.dp4 = function() {
  fl = list.files(vcf.directory, pattern = "lofreq-called-variants.vcf$", recursive = TRUE)
  out = vector("list", length(fl))
  for (i in seq_along(fl)) {
    samp = gsub("/.*", "", fl[i])
    ln = readLines(file.path(vcf.directory, fl[i]), warn = FALSE)
    ln = ln[!startsWith(ln, "#") & nzchar(ln)]
    if (!length(ln)) next
    f = data.table::tstrsplit(ln, "\t", fixed = TRUE)
    info = f[[8]]
    dp4 = sub(".*DP4=([0-9]+,[0-9]+,[0-9]+,[0-9]+).*", "\\1", info)
    dp4[!grepl("DP4=", info)] = NA_character_
    out[[i]] = data.table::data.table(sample = samp, locus = as.character(f[[1]]),
                                      position = as.numeric(f[[2]]), dp4_lo = dp4)
  }
  d = data.table::rbindlist(out)
  # FluLens keys the LoFreq record by position only and lets a later line win.
  if (nrow(d)) d = unique(d, by = c("sample", "locus", "position"), fromLast = TRUE)
  d
}
read.ivar.dp4 = function() {
  fl = list.files(vcf.directory, pattern = "ivar-called-variants.tsv$", recursive = TRUE)
  out = vector("list", length(fl))
  need = c("REGION", "POS", "ALT", "REF_DP", "REF_RV", "ALT_DP", "ALT_RV")
  for (i in seq_along(fl)) {
    samp = gsub("/.*", "", fl[i])
    d = try(data.table::fread(file.path(vcf.directory, fl[i]), sep = "\t", fill = TRUE,
                              colClasses = list(character = c("REGION", "REF", "ALT"))),
            silent = TRUE)
    if (inherits(d, "try-error") || is.null(nrow(d)) || nrow(d) == 0L ||
        !all(need %in% names(d))) next
    d = d[!grepl("^[+-]", ALT)]
    if (!nrow(d)) next
    out[[i]] = data.table::data.table(
      sample = samp, locus = as.character(d$REGION), position = as.numeric(d$POS),
      alternative = as.character(d$ALT),
      dp4_iv = paste(as.numeric(d$REF_DP) - as.numeric(d$REF_RV), as.numeric(d$REF_RV),
                     as.numeric(d$ALT_DP) - as.numeric(d$ALT_RV), as.numeric(d$ALT_RV),
                     sep = ","))
  }
  d = data.table::rbindlist(out)
  # iVar keys on the ALT base too, since it writes one row per alternative.
  if (nrow(d)) d = unique(d, by = c("sample", "locus", "position", "alternative"), fromLast = TRUE)
  d
}
lofreq.dp4 = read.lofreq.dp4()
ivar.dp4   = read.ivar.dp4()

fd = data.table::as.data.table(final.data)
orig.cols = names(final.data)   # merge() reorders; restored before the write
fd[, .ord := .I]
# Left-merge each DP4 source, then pick per caller. Both sources are unique on their
# key, so no row is duplicated. .ord restores the row order after the merges sort.
if (nrow(ivar.dp4)) {
  fd = merge(fd, ivar.dp4[, c("sample", "locus", "position", "alternative", "dp4_iv")],
             by = c("sample", "locus", "position", "alternative"), all.x = TRUE, sort = FALSE)
} else fd[, dp4_iv := NA_character_]
if (nrow(lofreq.dp4)) {
  fd = merge(fd, lofreq.dp4[, c("sample", "locus", "position", "dp4_lo")],
             by = c("sample", "locus", "position"), all.x = TRUE, sort = FALSE)
} else fd[, dp4_lo := NA_character_]
# iVar rows read iVar's TSV (keyed by alt); LoFreq and GATK4 rows read LoFreq's
# VCF at the position (GATK4 borrows it — same as FluLens strandRecOf).
fd[, dp4 := ifelse(method == "iVar", dp4_iv, dp4_lo)]
fd[, c("dp4_iv", "dp4_lo") := NULL]

d4 = data.table::tstrsplit(fd$dp4, ",", fixed = TRUE)
rf = as.numeric(d4[[1]]); rr = as.numeric(d4[[2]])
af = as.numeric(d4[[3]]); ar = as.numeric(d4[[4]])
at = af + ar; rt = rf + rr
alt.frac = ifelse(at > 0, af / at, 0)
ref.frac = ifelse(rt > 0, rf / rt, 0)
skew = abs(alt.frac - ref.frac)
freqF  = ifelse(is.finite(fd$allele_fraction), fd$allele_fraction, fd$allele_frequency)
depthN = suppressWarnings(as.numeric(fd$depth))
have   = !is.na(fd$dp4)

# strand class — the test order matters and matches assessCore.
strand = rep("not-assessed", nrow(fd))
if (run.strand.bias) {
  strand = ifelse(at < AS_STRAND_MIN_ALT, "too-few-alt",
           ifelse(alt.frac == 0 | alt.frac == 1, "skewed",
           ifelse(rt < AS_STRAND_MIN_ALT, "no-ref-control",
           ifelse(skew > AS_SKEW_BAD, "skewed",
           ifelse(skew > AS_SKEW_WARN, "some-skew", "balanced")))))
}
strand[!have] = NA_character_

# verdict codes: 0 real, 1 caution, 2 artefact, 3 cannot assess
verdict = rep(3L, nrow(fd))
is.bad  = have & run.strand.bias & strand == "skewed"
few.alt = have & at < AS_MIN_ALT
verdict[is.bad] = 2L
verdict[!is.bad & few.alt] = 2L
rem = have & verdict != 2L
caution = rem & (depthN < AS_MIN_DEPTH | freqF < AS_MIN_FREQ |
                 (run.strand.bias & strand %in% c("some-skew", "too-few-alt")))
verdict[caution] = 1L
verdict[rem & !caution] = 0L

fd$alt_reads    = at
fd$strand_class = strand
fd$assessment   = c("Looks real", "Treat with caution",
                    "Likely artefact", "Cannot assess")[verdict + 1L]

fd[, dp4 := NULL]   # internal only — its commas would break this unquoted CSV
data.table::setorder(fd, .ord)
fd[, .ord := NULL]
# Restore the original column order (merge moved the keys to the front); new columns
# go last, so the published layout does not change for a reader that uses position.
data.table::setcolorder(fd, c(orig.cols, "alt_reads", "strand_class", "assessment"))
final.data = as.data.frame(fd, stringsAsFactors = FALSE)

vt = table(factor(final.data$assessment,
                  levels = c("Looks real", "Treat with caution",
                             "Likely artefact", "Cannot assess")))
cat(sprintf("Call assessment (MIN_DEPTH=%g MIN_ALT=%g MIN_FREQ=%g strand_bias=%s): %s\n",
            AS_MIN_DEPTH, AS_MIN_ALT, AS_MIN_FREQ, run.strand.bias,
            paste(names(vt), vt, sep = "=", collapse = "  ")))

#Saves the data
write.csv(final.data, paste0(output.directory, "/", save.name, ".csv"),
          row.names = F, quote = F)




#########################
###### END SCRIPT
#########################














