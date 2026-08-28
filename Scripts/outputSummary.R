#### Prior to this script, run the Flumina pipeline to obtain variant calls in VCF files
#### Next, run 2_convertVCFtoTable.R to get an amino acid table, all_sample_amino_acids.txt

#Takes seconds to run on 500 samples

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

#Define these

# curated csv database with the columns "Gene", "Amino_Acid", and "Type" (category of site) of interest
# Other columns can be added and joined with the variants
aa.table.path = gsub("\"", "", config$AA_DB)

# The curated database is optional. Everything this script produces is a join
# against it, so with no database there is nothing to summarise — but the full
# variant table and the complete amino-acid table have already been written by
# the earlier steps, and those are the substantive outputs. Stop cleanly rather
# than failing the run over a file the user chose not to supply.
if (length(aa.table.path) == 0L || aa.table.path == "" ||
    aa.table.path == "NULL" || !file.exists(aa.table.path)) {
  message("No curated amino acid database provided (AA_DB); ",
          "skipping curated_amino_acids.txt and summary_curated_sites.txt.")
  message("The full variant table and all_sample_amino_acids.txt are unaffected.")
  quit(save = "no", status = 0)
}

# output.directory used from 1_convertVCFtoTable.R
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/variant_analysis")

#Grouping category column name joined in step 2 with the variant sample data
#Set to NULL if there are no groupings to use
group.names = gsub("\"", "", config$GROUP_NAMES)
if(length(group.names) == 0L || group.names == "NULL") {
  group.names <- NULL
}

#gets filters
depth.val = gsub("\"", "", config$MIN_DEPTH)
qual.val = gsub("\"", "", config$MIN_QUALITY)
af.val = gsub("\"", "", config$MIN_ALLELE_FREQUENCY)


#### debuggin
#output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Bird_Flu/bird_flu_new/variant_analysis"
#aa.table.path = paste0("/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Bird_Flu/curated_database.csv")
#threads = 4
#group.names = "discrete_host"   # a COLUMN of the METADATA csv, not a path

#############################################
#### Should not need to modify below here
#############################################

#output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Bird_Flu/variant_analysis"

#read in previous database
sample.data = read.table(paste0(output.directory, "/all_sample_amino_acids.txt"), sep = "\t", header = T, na.strings = "")

#GROUP_NAMES is the name of a metadata COLUMN (e.g. "discrete_host"), not a path.
#findAAChanges.R already merged the METADATA csv into all_sample_amino_acids.txt,
#so the column is present here and no second file needs to be read. The grouping
#loop below keys off sample.data$group, so map the chosen column onto it.
if (is.null(group.names) != TRUE){
  if (!group.names %in% colnames(sample.data)){
    stop(paste0("GROUP_NAMES column '", group.names, "' not found in ",
                output.directory, "/all_sample_amino_acids.txt.\n",
                "  Available columns: ", paste(colnames(sample.data), collapse = ", "), "\n",
                "  GROUP_NAMES must name a column of the METADATA csv, or be NULL."))
  }
  sample.data$group = sample.data[[group.names]]
} #end if

#Apply the config depth/quality floor (MIN_DEPTH / MIN_QUALITY). The amino-acid
#table from findAAChanges.R is already filtered; this re-applies the same
#thresholds defensively so summaries always honor the floor even if run on an
#older table. Allele-frequency handling is left to the summary section below.
if (length(depth.val) != 0L && nzchar(depth.val)){
  sample.data = sample.data[!is.na(sample.data$depth) & sample.data$depth >= as.numeric(depth.val),]
}
if (length(qual.val) != 0L && nzchar(qual.val)){
  sample.data = sample.data[!is.na(sample.data$quality) & sample.data$quality >= as.numeric(qual.val),]
}


# if (group.names == "AUTO"){
#   
#   sample.names = unique(sample.data$sample)
#   full.name = gsub("^[^_]*_", "", sample.names)
#   full.name = gsub("_.*", "", full.name)
#   full.name[full.name %in% names(table(full.name)[table(full.name) <= 2])] = "Wild-Bird"
#   name.data = data.frame(sample = sample.names, discrete_host = full.name)
#   
#   #Save large tab delimited table of all the amino acids
#   write.csv(name.data, paste0(output.directory, "/sample_names.csv"), quote = F)
#   
# }

#Reads in curated database
best.aa = read.csv(aa.table.path, header = TRUE, sep = ",")
# A UTF-8 BOM in the curated CSV mangles the first column name (e.g. "X...Gene"),
# which silently nulls best.aa$Gene and makes the merge below return 0 rows.
# Normalise the first column back to "Gene".
colnames(best.aa)[1] = "Gene"
best.aa[is.na(best.aa) == TRUE] = "NA"

#############################################
#### Which reading frame the curated positions are in
#############################################
# The amino-acid table now carries one row per PRODUCT, so a single segment can
# contribute rows in two reading frames (PA and PA-X, M1 and M2, NS1 and NEP).
# The curated database keys on the SEGMENT name, and its coordinates are
# primary-product coordinates: it lists M1 N30D and M1 T215A as virulence
# determinants, and NS1 P42S and D92E, all of which are primary-ORF residues.
#
# Joining it against secondary products as well would add 22 "hits" on the
# swine WGS data (MP 30/43/77 onto M2, NS 42/91/92/94 onto NEP). Every one is a
# numeric coincidence, not a marker — so restrict the join to primary products.
#
# A curated database that genuinely wants a secondary-product site can say so
# by adding a `Product` column; when present it is matched against `product`
# instead, and segment-keyed rows keep their primary-only behaviour.
if ("product_primary" %in% colnames(sample.data)){
  has.product.col = "Product" %in% colnames(best.aa)
  if (!has.product.col){
    n.before = nrow(sample.data)
    sample.data = sample.data[!is.na(sample.data$product_primary) &
                              sample.data$product_primary %in% c(TRUE, "TRUE"), ]
    message(sprintf("Curated join restricted to primary products: %d -> %d rows. %s",
                    n.before, nrow(sample.data),
                    "Add a 'Product' column to the curated CSV to target M2 / NEP / PA-X / PB1-F2."))
  }
}

#Obtains gene names
gene.names = unique(sample.data$locus)

#Combines all the data into 1 data frame
save.sample = c()
for (i in seq_along(gene.names)){

  #Subsets data
  temp.sample = sample.data[sample.data$locus %in% gene.names[i],]
  temp.fun = best.aa[best.aa$Gene %in% gene.names[i],]

  if (nrow(temp.fun) == 0){ next }

  #Merges with curated database. With a Product column the frame is explicit,
  #so match on it as well as the position.
  if ("Product" %in% colnames(best.aa) && "product" %in% colnames(temp.sample)){
    merge.sample = merge(temp.sample, temp.fun,
                         by.x = c("product", "aa_position"),
                         by.y = c("Product",  "Amino_Acid"))
  } else {
    merge.sample = merge(temp.sample, temp.fun, by.x = "aa_position", by.y = "Amino_Acid")
  }
  save.sample = rbind(save.sample, merge.sample)

} #end i

#Save large tab delimited table of all the amino acids
write.table(save.sample, paste0(output.directory, "/curated_amino_acids.txt"),
            row.names = F, quote = F, sep = "\t")


################################
### Make summary table
################################

#Only amino acid changing sites
red.samples = save.sample[save.sample$aa_changing == "YES",]

#Removes duplicate stuff so they aren't counted twice
red.samples$sample = gsub("-r_", "_", red.samples$sample)
red.samples$sample = gsub("-v_", "_", red.samples$sample)

#minimum allele frequency 
red.samples = red.samples[red.samples$allele_frequency >= 0.005,]

#Obtains different animal groups, if desired
if (is.null(group.names) != TRUE){
  group.values = unique(red.samples$group)  
} else{ group.values = "all" }

all.data = c()
for (i in seq_along(group.values)){

  #Obtains different animal groups
  if (is.null(group.names) != TRUE){
    group.data = red.samples[red.samples$group %in% group.values[i],]
  } else{ group.data = red.samples }  
  
  #Obtains gene names
  gene.names = unique(group.data$locus)
  
  #Creates table
  sample.table = c()
  for (j in seq_along(gene.names)){
    
    #Subsets data
    gene.data = group.data[grep(gene.names[j], group.data$locus),]
    
    #Gets unique positions
    aa.pos = unique(gene.data$aa_position)
    
    #empty table
    temp.table = data.frame(animal = as.character(),
                            gene = as.character(),
                            aa_position = as.numeric(),
                            mutation = as.character(),
                            funct = as.character(),
                            no_animal = as.numeric(),
                            ave_alle_freq = as.numeric(),
                            consensus = as.numeric(),
                            low_freq = as.numeric(),
                            gatk4 = as.numeric())
    
    #Loops through positions
    for (k in seq_along(aa.pos)){
      
      #empty data for this section
      vector.table = data.frame(animal = as.character(),
                                gene = as.character(),
                                aa_position = as.numeric(),
                                mutation = as.character(),
                                funct = as.character(),
                                no_animal = as.numeric(),
                                ave_alle_freq = as.numeric(),
                                consensus = as.numeric(),
                                low_freq = as.numeric(),
                                gatk4 = as.numeric())
      
      #Subsets to amino acid data
      aa.data = gene.data[gene.data$aa_position == aa.pos[k],]
      lf.data = aa.data[aa.data$method %in% "LoFreq",]
      gk.data = aa.data[aa.data$method %in% "GATK4",]
      
      ##### Groupings categories
      vector.table[1,1] = group.values[i]
      vector.table[1,2] = gene.names[j]
      vector.table[1,3] = aa.pos[k]
      
      #Formats mutation names
      ref.chars = unique(aa.data$reference_aa)
      aa.chars = unique(aa.data$alternative_aa)
      new.chars = paste0(ref.chars, aa.pos[k], aa.chars)
      vector.table[1,4] = paste0(new.chars, collapse = ", ")
      
      #saves AA type
      vector.table[1,5] = unique(aa.data$Type)
      
      #Counts samples
      vector.table[1,6] = length(unique(aa.data$sample))
      
      #mean allele frequency
      vector.table[1,7] = round(mean(lf.data$allele_frequency), 3)
      
      #gets those in consensus
      freq.aa = aa.data[aa.data$allele_frequency > 0.50,]
      vector.table[1,8] = length(unique(freq.aa$sample))
      
      #gets those in low frequency
      freq.aa = lf.data[lf.data$allele_frequency < 0.50,]
      vector.table[1,9] = length(unique(freq.aa$sample))
      
      #gets those in gatk4
      vector.table[1,10] = length(unique(gk.data$sample))
      
      #combines togther
      temp.table = rbind(temp.table, vector.table)
      
    }# end k loop
    
    #saves sample data
    sample.table = rbind(sample.table, temp.table)
    
  }#end j loop
    
  #saves all the data
  all.data = rbind(all.data, sample.table)
  
}# i loop

#Saves final data
write.table(all.data, paste0(output.directory, "/summary_curated_sites.txt"),
            row.names = F, quote = F, sep = "\t")












