#### Run the Flumina pipeline first to obtain variant calls in VCF files.
#### Then run convertVCFtoTable.R to create all_sample_amino_acids.txt.

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

# Define these parameters.

# curated csv database with the columns "Gene", "Amino_Acid", and "Type" (category of site) of interest
# Other columns can be added and joined with the variants
aa.table.path = gsub("\"", "", config$AA_DB)

# The curated database is optional. This script produces joins against it, so
# there is no curated summary without the database. The full variant and
# amino-acid tables were already written by earlier steps; stop cleanly rather
# than failing because an optional file was omitted.
if (length(aa.table.path) == 0L || aa.table.path == "" ||
    aa.table.path == "NULL" || !file.exists(aa.table.path)) {
  message("No curated amino acid database provided (AA_DB); ",
          "skipping curated_amino_acids.txt and summary_curated_sites.txt.")
  message("The full variant table and all_sample_amino_acids.txt are unaffected.")
  quit(save = "no", status = 0)
}

# output.directory used from convertVCFtoTable.R
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/variant_analysis")

# Grouping column joined to the variant sample data in step 2.
# Set to NULL when no grouping is required.
group.names = gsub("\"", "", config$GROUP_NAMES)
if(length(group.names) == 0L || group.names == "NULL") {
  group.names <- NULL
}

# Read the configured filters.
depth.val = gsub("\"", "", config$MIN_DEPTH)
qual.val = gsub("\"", "", config$MIN_QUALITY)
af.val = gsub("\"", "", config$MIN_ALLELE_FREQUENCY)


#### Debugging examples
#output.directory = "path/to/variant_analysis"
#aa.table.path = paste0("path/to/curated_database.csv")
#threads = 4
#group.names = "discrete_host"   # a COLUMN of the METADATA csv, not a path

#############################################
#### Should not need to modify below here
#############################################

# Read the previously generated amino-acid table.
sample.data = read.table(paste0(output.directory, "/all_sample_amino_acids.txt"), sep = "\t", header = T, na.strings = "")

# GROUP_NAMES is the name of a metadata column (for example, "discrete_host"),
# not a path. findAAChanges.R already merged the metadata CSV into
# all_sample_amino_acids.txt, so no second file is required. The grouping loop
# uses sample.data$group, so map the configured column to it.
if (is.null(group.names) != TRUE){
  if (!group.names %in% colnames(sample.data)){
    stop(paste0("GROUP_NAMES column '", group.names, "' not found in ",
                output.directory, "/all_sample_amino_acids.txt.\n",
                "  Available columns: ", paste(colnames(sample.data), collapse = ", "), "\n",
                "  GROUP_NAMES must name a column of the METADATA csv, or be NULL."))
  }
  sample.data$group = sample.data[[group.names]]
} # End of if block.

#Apply the config depth/quality floor (MIN_DEPTH / MIN_QUALITY). The amino-acid
# table from findAAChanges.R is already filtered; this reapplies the same
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
#   # Save the complete tab-delimited amino-acid table.
#   write.csv(name.data, paste0(output.directory, "/sample_names.csv"), quote = F)
#   
# }

# Read the curated database.
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
# Joining it against secondary products as well would add spurious hits — a
# curated position matching by number, not by biology — so restrict the join to
# primary products.
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

# Obtain gene names.
gene.names = unique(sample.data$locus)

# Combine all records into a single data frame.
save.sample = c()
for (i in seq_along(gene.names)){

  # Subset the data.
  temp.sample = sample.data[sample.data$locus %in% gene.names[i],]
  temp.fun = best.aa[best.aa$Gene %in% gene.names[i],]

  if (nrow(temp.fun) == 0){ next }

  # Merge with the curated database. When a Product column is present, the
  # reading frame is explicit, so match it as well as the position.
  if ("Product" %in% colnames(best.aa) && "product" %in% colnames(temp.sample)){
    merge.sample = merge(temp.sample, temp.fun,
                         by.x = c("product", "aa_position"),
                         by.y = c("Product",  "Amino_Acid"))
  } else {
    merge.sample = merge(temp.sample, temp.fun, by.x = "aa_position", by.y = "Amino_Acid")
  }
  save.sample = rbind(save.sample, merge.sample)

} # End of i loop.

# Save the complete tab-delimited amino-acid table.
write.table(save.sample, paste0(output.directory, "/curated_amino_acids.txt"),
            row.names = F, quote = F, sep = "\t")


################################
### Make summary table
################################

# Retain only amino-acid-changing sites.
red.samples = save.sample[save.sample$aa_changing == "YES",]

# Remove duplicate records so they are not counted twice.
red.samples$sample = gsub("-r_", "_", red.samples$sample)
red.samples$sample = gsub("-v_", "_", red.samples$sample)

# Apply the minimum allele-frequency threshold.
red.samples = red.samples[red.samples$allele_frequency >= 0.005,]

# Obtain the requested animal groups, when configured.
if (is.null(group.names) != TRUE){
  group.values = unique(red.samples$group)  
} else{ group.values = "all" }

all.data = c()
for (i in seq_along(group.values)){

  # Obtain the distinct groups.
  if (is.null(group.names) != TRUE){
    group.data = red.samples[red.samples$group %in% group.values[i],]
  } else{ group.data = red.samples }  
  
  # Obtain gene names.
  gene.names = unique(group.data$locus)
  
  # Initialize the summary table.
  sample.table = c()
  for (j in seq_along(gene.names)){
    
    # Subset the data.
    gene.data = group.data[grep(gene.names[j], group.data$locus),]
    
    # Obtain unique positions.
    aa.pos = unique(gene.data$aa_position)
    
    # Initialize an empty table.
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
    
    # Iterate over positions.
    for (k in seq_along(aa.pos)){
      
      # Initialize the data for this section.
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
      
      # Subset the amino-acid data.
      aa.data = gene.data[gene.data$aa_position == aa.pos[k],]
      lf.data = aa.data[aa.data$method %in% "LoFreq",]
      gk.data = aa.data[aa.data$method %in% "GATK4",]
      
      ##### Grouping categories
      vector.table[1,1] = group.values[i]
      vector.table[1,2] = gene.names[j]
      vector.table[1,3] = aa.pos[k]
      
      # Format mutation names.
      ref.chars = unique(aa.data$reference_aa)
      aa.chars = unique(aa.data$alternative_aa)
      new.chars = paste0(ref.chars, aa.pos[k], aa.chars)
      vector.table[1,4] = paste0(new.chars, collapse = ", ")
      
      # Store the amino-acid type.
      vector.table[1,5] = unique(aa.data$Type)
      
      # Count samples.
      vector.table[1,6] = length(unique(aa.data$sample))
      
      # Calculate the mean allele frequency.
      vector.table[1,7] = round(mean(lf.data$allele_frequency), 3)
      
      # Identify consensus calls.
      freq.aa = aa.data[aa.data$allele_frequency > 0.50,]
      vector.table[1,8] = length(unique(freq.aa$sample))
      
      # Identify low-frequency calls.
      freq.aa = lf.data[lf.data$allele_frequency < 0.50,]
      vector.table[1,9] = length(unique(freq.aa$sample))
      
      # Identify GATK4 calls.
      vector.table[1,10] = length(unique(gk.data$sample))
      
      # Combine the summary fields.
      temp.table = rbind(temp.table, vector.table)
      
    }# end k loop
    
    # Store the sample-level data.
    sample.table = rbind(sample.table, temp.table)
    
  } # End of j loop.
    
  # Store the complete summary.
  all.data = rbind(all.data, sample.table)
  
}# i loop

# Save the final output.
write.table(all.data, paste0(output.directory, "/summary_curated_sites.txt"),
            row.names = F, quote = F, sep = "\t")







