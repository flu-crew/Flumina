#### makeGTF.R
#### Generates per-segment GTF annotation files from a reference FASTA for use
#### with SNPGenie.  Identifies influenza A segment types from sequence names,
#### verifies start/stop codon positions, and writes individual GTF files plus
#### per-segment FASTAs.  A combined GTF is also written when the reference
#### contains more than one segment.

args = commandArgs(trailingOnly = TRUE)
#args = "config.cfg"

# Parse configuration file
lines <- readLines(args)
config <- list()
for (line in lines) {
  line <- trimws(line)
  if (nchar(line) != 0 && !startsWith(line, "#")) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key   <- trimws(parts[1])
      value <- trimws(parts[2])
      config[[key]] <- value
    }
  }
}

reference.path   = gsub("\"", "", config$REFERENCE_FILE)
output.directory = paste0(gsub("\"", "", config$OUTPUT_DIRECTORY), "/reference_gtf")

dir.create(output.directory, showWarnings = FALSE)

script.dir = dirname(sub("^--file=", "",
                         grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(script.dir) || !nzchar(script.dir)) script.dir = "."
source(file.path(script.dir, "fluORFs.R"))

require(Biostrings)
ref.seqs = readDNAStringSet(reference.path)

#############################################
#### Helper functions
#############################################

# Identify the influenza A segment type from a sequence name.
# Delegates to fluORFs.R so the segment rules exist in exactly one place.
get_segment_type <- function(name) flu_segment_type(name)

# Standard product names for influenza A genes
gene_product <- function(gene) {
  products = c(
    HA = "hemagglutinin",         "NA" = "neuraminidase",
    NP = "nucleocapsid protein",  PB1 = "polymerase PB1",
    PB2 = "polymerase PB2",       PA  = "polymerase PA",
    "PB1-F2" = "PB1-F2 protein",  "PA-X" = "PA-X protein",
    M1  = "matrix protein 1",     M2  = "matrix protein 2",
    NS1 = "nonstructural protein 1", NEP = "nuclear export protein"
  )
  if (gene %in% names(products)) return(products[[gene]])
  return(gene)
}

# Assemble a single tab-delimited GTF row
gtf_row <- function(seqname, feature, start, end,
                    frame = ".", strand = "+", attributes) {
  paste(seqname, ".", feature, start, end, ".", strand, frame,
        attributes, sep = "\t")
}

# Build a GTF attribute string from named key=value pairs
attr_str <- function(...) {
  pairs = list(...)
  parts = mapply(function(k, v) paste0(k, ' "', v, '"'),
                 names(pairs), unlist(pairs))
  paste(parts, collapse = "; ")
}

# Verify a codon in the sequence at a given 1-based position
check_codon <- function(seq.str, pos, type = "start", label = "") {
  if (pos < 1 || pos + 2 > nchar(seq.str)) {
    warning(paste0("  [", label, "] codon position ", pos, " out of range"))
    return(FALSE)
  }
  codon = substr(seq.str, pos, pos + 2)
  if (type == "start" && codon != "ATG") {
    warning(paste0("  [", label, "] Expected ATG start codon at pos ", pos,
                   " but found: ", codon))
    return(FALSE)
  }
  if (type == "stop" && !codon %in% c("TAA", "TAG", "TGA")) {
    warning(paste0("  [", label, "] Expected stop codon at pos ", pos,
                   " but found: ", codon))
    return(FALSE)
  }
  return(TRUE)
}

#############################################
#### GTF entry builders
####
#### The CDS intervals come from Scripts/fluORFs.R — the same function the
#### variant pipeline uses to turn a nucleotide position into an amino-acid
#### position. They used to be written out twice, here and there, with a comment
#### asking whoever edited one to remember the other. They diverged: fluORFs
#### gained stop-codon trimming (NS1 and PA-X have strain-variable C-termini —
#### 219 vs 230 aa and 232 vs 252 aa between the H3N2 and H5N1 references) and
#### this file did not. One source now, so it cannot happen again.
#############################################

flu_gene_biotype <- "protein_coding"

make_gtf_entries <- function(seqname, seq.len, seg.type, seq.str) {

  if (is.na(seg.type)) {
    warning(paste0("Could not identify segment type for '", seqname,
                   "' — writing single-CDS annotation."))
  }

  orfs = flu_orfs_for(seqname, seq.str, seg.type)
  # NCBI lists the spliced product first on MP and NS; keep that convention so
  # the emitted GTF still matches a reference annotation line for line.
  if (!is.na(seg.type) && seg.type %in% c("MP", "NS")) orfs = rev(orfs)
  rows = c()
  cat("  Segment:", seqname, "| genes:",
      paste(sapply(orfs, `[[`, "gene"), collapse = " + "), "| len:", seq.len, "\n")

  for (o in orfs) {
    gene   = o$gene
    ex     = o$exons
    n.ex   = nrow(ex)
    cds.end   = ex$end[n.ex]
    stop.start= cds.end + 1L
    stop.end  = min(cds.end + 3L, seq.len)
    tx = paste0("unassigned_transcript_", which(sapply(orfs, `[[`, "gene") == gene)[1])

    check_codon(seq.str, ex$start[1], "start", gene)
    if (stop.end - stop.start == 2L) check_codon(seq.str, stop.start, "stop", gene)

    # Per-gene annotation NCBI carries and SNPGenie ignores, but which makes the
    # output comparable to a downloaded reference GTF.
    # Two slots, because NCBI puts them in different places: PA-X's `exception`
    # sits before gbkey, NEP's `note` after gene.
    extra = if (gene == "PA-X") list(exception = "ribosomal slippage") else list()
    mid   = if (gene == "NEP")  list(note = "nonstructural protein 2") else list()

    gene.attrs = c(list(gene_id = gene, transcript_id = "", gbkey = "Gene",
                        gene = gene, gene_biotype = flu_gene_biotype),
                   if (gene == "NEP") list(gene_synonym = "NS2") else list())
    rows = c(rows,
      gtf_row(seqname, "gene", ex$start[1], stop.end, ".", "+",
        do.call(attr_str, gene.attrs)))

    before = 0L
    for (k in seq_len(n.ex)) {
      # GTF frame: bases to skip at the start of this exon to reach a codon
      frame = (3L - (before %% 3L)) %% 3L
      a = do.call(attr_str, c(list(gene_id = gene, transcript_id = tx), extra,
                   list(gbkey = "CDS", gene = gene), mid,
                   if (n.ex > 1) list(part = as.character(k)) else list(),
                   list(product = gene_product(gene), exon_number = as.character(k))))
      rows = c(rows, gtf_row(seqname, "CDS", ex$start[k], ex$end[k],
                             as.character(frame), "+", a))
      before = before + (ex$end[k] - ex$start[k] + 1L)
    }

    rows = c(rows,
      gtf_row(seqname, "start_codon", ex$start[1], ex$start[1] + 2L, "0", "+",
        do.call(attr_str, c(list(gene_id = gene, transcript_id = tx), extra,
                 list(gbkey = "CDS", gene = gene), mid,
                 list(product = gene_product(gene), exon_number = "1")))),
      gtf_row(seqname, "stop_codon", stop.start, stop.end, "0", "+",
        do.call(attr_str, c(list(gene_id = gene, transcript_id = tx), extra,
                 list(gbkey = "CDS", gene = gene), mid,
                 list(product = gene_product(gene),
                      exon_number = as.character(n.ex))))))
  }
  rows
}

# Re-offset all coordinate columns in a vector of GTF entry strings and
# optionally rename the chromosome (first field) to a new name.
# Header/terminator lines (#gtf-version, ###) are passed through unchanged.
offset_gtf_entries <- function(entries, offset, new_chrom = NULL) {
  sapply(entries, function(line) {
    if (grepl("^#", line)) return(line)
    fields        = strsplit(line, "\t")[[1]]
    if (!is.null(new_chrom)) fields[1] = new_chrom
    fields[4]     = as.character(as.integer(fields[4]) + offset)  # start
    fields[5]     = as.character(as.integer(fields[5]) + offset)  # end
    paste(fields, collapse = "\t")
  }, USE.NAMES = FALSE)
}

#############################################
#### Main: process each sequence in the reference FASTA
#############################################

cat("Creating GTF files in:", output.directory, "\n")

combined.entries  = c()
combined.seq.str  = ""
cumulative.offset = 0L
offsets.data      = data.frame(segment = character(), offset = integer(),
                               length  = integer(),   stringsAsFactors = FALSE)

for (i in seq_along(ref.seqs)) {

  seq.name = names(ref.seqs)[i]
  seq.str  = as.character(ref.seqs[[i]])
  seq.len  = nchar(seq.str)
  seg.type = get_segment_type(seq.name)

  entries = make_gtf_entries(seq.name, seq.len, seg.type, seq.str)

  # Accumulate offset entries and sequence for the combined files
  offset.entries    = offset_gtf_entries(entries, cumulative.offset, new_chrom = "combined")
  combined.entries  = c(combined.entries, offset.entries)
  combined.seq.str  = paste0(combined.seq.str, seq.str)
  offsets.data      = rbind(offsets.data,
                            data.frame(segment = seq.name,
                                       offset  = cumulative.offset,
                                       length  = seq.len,
                                       stringsAsFactors = FALSE))
  cumulative.offset = cumulative.offset + seq.len

  # Write per-segment GTF
  out.gtf = paste0(output.directory, "/", seq.name, ".gtf")
  writeLines(c("#gtf-version 2.2", entries, "###"), out.gtf)
  cat("  Wrote:", out.gtf, "\n")

  # Write per-segment FASTA (required by SNPGenie alongside the GTF)
  out.fasta = paste0(output.directory, "/", seq.name, ".fasta")
  writeXStringSet(ref.seqs[i], out.fasta)
  cat("  Wrote:", out.fasta, "\n")
}

# Write combined files only when the reference contains multiple segments
if (length(ref.seqs) > 1) {

  # Combined GTF with offset coordinates
  out.combined.gtf = paste0(output.directory, "/combined.gtf")
  writeLines(c("#gtf-version 2.2", combined.entries, "###"), out.combined.gtf)
  cat("Wrote combined GTF:", out.combined.gtf, "\n")

  # Combined FASTA (all segments concatenated into one sequence named "combined")
  combined.dna = Biostrings::DNAStringSet(Biostrings::DNAString(combined.seq.str))
  names(combined.dna) = "combined"
  out.combined.fasta = paste0(output.directory, "/combined.fasta")
  writeXStringSet(combined.dna, out.combined.fasta)
  cat("Wrote combined FASTA:", out.combined.fasta, "\n")

  # Offsets table (segment order, offset, and length) used by runSNPGenie.R
  out.offsets = paste0(output.directory, "/combined_offsets.tsv")
  write.table(offsets.data, out.offsets,
              row.names = FALSE, quote = FALSE, sep = "\t")
  cat("Wrote offsets table:", out.offsets, "\n")
}
cat("GTF creation complete.\n")
