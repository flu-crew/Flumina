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

for (i in seq_along(ivar.files)) {

  ivar.path = paste0(vcf.directory, "/", ivar.files[i])
  ivar.tab = try(utils::read.table(ivar.path, sep = "\t", header = TRUE,
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
    allele_frequency = numeric(), aa_position = numeric(), ivar_pass = character())
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

# ivar_pass belongs only to iVar rows, so the other two callers are given the
# column as NA before the bind. NA means "this caller has no such verdict",
# which is a different statement from FALSE and must not collapse into it.
lofreq.data$ivar_pass = NA_character_
collect.data$ivar_pass = NA_character_

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

#Saves the data
write.csv(final.data, paste0(output.directory, "/", save.name, ".csv"),
          row.names = F, quote = F)




#########################
###### END SCRIPT
#########################














