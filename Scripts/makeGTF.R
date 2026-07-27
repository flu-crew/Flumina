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

require(Biostrings)
ref.seqs = readDNAStringSet(reference.path)

#############################################
#### Helper functions
#############################################

# Identify the influenza A segment type from a sequence name
get_segment_type <- function(name) {
  n = toupper(name)
  if (grepl("PB2", n))                            return("PB2")
  if (grepl("PB1", n))                            return("PB1")
  if (grepl("NP",  n) && !grepl("PB", n))         return("NP")
  if (grepl("PA",  n) && !grepl("PB", n))         return("PA")
  if (grepl("HA|HEMA|H[0-9]", n))                 return("HA")
  if (grepl("NA|N[0-9]", n) && !grepl("NP|PB|PAND", n)) return("NA")
  if (grepl("MP|MATRIX|_M$|_M[12_]", n))          return("MP")
  if (grepl("NS", n))                              return("NS")
  return(NA)
}

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

# Proportionally scale a canonical coordinate to a different sequence length
scale_pos <- function(pos, canonical_len, actual_len) {
  round(pos * actual_len / canonical_len)
}

#############################################
#### GTF entry builders per segment type
####
#### Canonical coordinates come from verified H3N2 influenza A references.
#### For single-gene segments the CDS boundaries are derived directly from
#### the actual sequence length.  For spliced / alternate-ORF segments the
#### internal coordinates are scaled proportionally when the input length
#### differs from the canonical reference.
####
#### NOTE: the CDS intervals below are also expressed, independently, in
#### Scripts/fluORFs.R, which is what convertVCFtoTable.R / findAAChanges.R /
#### runWFABC.R use to turn a nucleotide position into an amino-acid position.
#### The two must agree. They are kept separate because this script also emits
#### gene / start_codon / stop_codon rows that the position mapping has no use
#### for — but if you change a coordinate here, change it there too, and run
#### Scripts/test_fluORFs.R, which fails when they diverge.
#############################################

make_gtf_entries <- function(seqname, seq.len, seg.type, seq.str) {

  rows = c()

  if (is.na(seg.type)) {
    warning(paste0("Could not identify segment type for '", seqname,
                   "' — writing single-CDS annotation."))
    seg.type = "UNKNOWN"
  }

  # ------------------------------------------------------------------
  # Single-gene segments: HA, NA, NP, PB2, UNKNOWN
  # CDS is always 1 to (len-3); stop codon is last 3 bases.
  # ------------------------------------------------------------------
  if (seg.type %in% c("HA", "NA", "NP", "PB2", "UNKNOWN")) {

    gene.name  = if (seg.type == "UNKNOWN") seqname else seg.type
    cds.end    = seq.len - 3
    stop.start = seq.len - 2

    cat("  Segment:", seqname, "| gene:", gene.name, "| len:", seq.len, "\n")
    check_codon(seq.str, 1,          "start", gene.name)
    check_codon(seq.str, stop.start, "stop",  gene.name)

    rows = c(rows,
      gtf_row(seqname, "gene", 1, seq.len, ".", "+",
        attr_str(gene_id = gene.name, transcript_id = "",
                 gbkey = "Gene", gene = gene.name, gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, cds.end, "0", "+",
        attr_str(gene_id = gene.name, transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = gene.name,
                 product = gene_product(gene.name), exon_number = "1")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = gene.name, transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = gene.name,
                 product = gene_product(gene.name), exon_number = "1")),
      gtf_row(seqname, "stop_codon", stop.start, seq.len, "0", "+",
        attr_str(gene_id = gene.name, transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = gene.name,
                 product = gene_product(gene.name), exon_number = "1"))
    )
  }

  # ------------------------------------------------------------------
  # PB1: primary ORF + PB1-F2 alternate ORF (canonical len 2274)
  # PB1-F2 starts in a different reading frame at position 95.
  # ------------------------------------------------------------------
  else if (seg.type == "PB1") {

    canonical_len  = 2274
    f2.start       = scale_pos(95,  canonical_len, seq.len)
    f2.cds.end     = scale_pos(364, canonical_len, seq.len)
    f2.stop.start  = scale_pos(365, canonical_len, seq.len)
    f2.stop.end    = scale_pos(367, canonical_len, seq.len)

    cat("  Segment:", seqname, "| genes: PB1 + PB1-F2 | len:", seq.len, "\n")
    if (seq.len != canonical_len)
      cat("  Note: scaling PB1-F2 coordinates from canonical",
          canonical_len, "to", seq.len, "\n")
    check_codon(seq.str, 1,             "start", "PB1")
    check_codon(seq.str, seq.len - 2,   "stop",  "PB1")
    check_codon(seq.str, f2.start,      "start", "PB1-F2")
    check_codon(seq.str, f2.stop.start, "stop",  "PB1-F2")

    rows = c(rows,
      # PB1 primary ORF
      gtf_row(seqname, "gene", 1, seq.len, ".", "+",
        attr_str(gene_id = "PB1", transcript_id = "",
                 gbkey = "Gene", gene = "PB1", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, seq.len - 3, "0", "+",
        attr_str(gene_id = "PB1", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PB1",
                 product = gene_product("PB1"), exon_number = "1")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "PB1", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PB1",
                 product = gene_product("PB1"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", seq.len - 2, seq.len, "0", "+",
        attr_str(gene_id = "PB1", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PB1",
                 product = gene_product("PB1"), exon_number = "1")),
      # PB1-F2 alternate ORF
      gtf_row(seqname, "gene", f2.start, f2.stop.end, ".", "+",
        attr_str(gene_id = "PB1-F2", transcript_id = "",
                 gbkey = "Gene", gene = "PB1-F2", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", f2.start, f2.cds.end, "0", "+",
        attr_str(gene_id = "PB1-F2", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "PB1-F2",
                 product = gene_product("PB1-F2"), exon_number = "1")),
      gtf_row(seqname, "start_codon", f2.start, f2.start + 2, "0", "+",
        attr_str(gene_id = "PB1-F2", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "PB1-F2",
                 product = gene_product("PB1-F2"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", f2.stop.start, f2.stop.end, "0", "+",
        attr_str(gene_id = "PB1-F2", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "PB1-F2",
                 product = gene_product("PB1-F2"), exon_number = "1"))
    )
  }

  # ------------------------------------------------------------------
  # PA: primary ORF + PA-X ribosomal-slippage product (canonical len 2151)
  # PA-X exon1 shares sequence with PA; exon2 is in the +1 reading frame.
  # ------------------------------------------------------------------
  else if (seg.type == "PA") {

    canonical_len   = 2151
    px.e1.end       = scale_pos(570, canonical_len, seq.len)
    px.e2.start     = scale_pos(572, canonical_len, seq.len)
    px.e2.end       = scale_pos(697, canonical_len, seq.len)
    px.stop.start   = scale_pos(698, canonical_len, seq.len)
    px.stop.end     = scale_pos(700, canonical_len, seq.len)

    cat("  Segment:", seqname, "| genes: PA + PA-X | len:", seq.len, "\n")
    if (seq.len != canonical_len)
      cat("  Note: scaling PA-X coordinates from canonical",
          canonical_len, "to", seq.len, "\n")
    check_codon(seq.str, 1,           "start", "PA")
    check_codon(seq.str, seq.len - 2, "stop",  "PA")

    rows = c(rows,
      # PA primary ORF
      gtf_row(seqname, "gene", 1, seq.len, ".", "+",
        attr_str(gene_id = "PA", transcript_id = "",
                 gbkey = "Gene", gene = "PA", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, seq.len - 3, "0", "+",
        attr_str(gene_id = "PA", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PA",
                 product = gene_product("PA"), exon_number = "1")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "PA", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PA",
                 product = gene_product("PA"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", seq.len - 2, seq.len, "0", "+",
        attr_str(gene_id = "PA", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "PA",
                 product = gene_product("PA"), exon_number = "1")),
      # PA-X ribosomal slippage product
      gtf_row(seqname, "gene", 1, px.stop.end, ".", "+",
        attr_str(gene_id = "PA-X", transcript_id = "",
                 gbkey = "Gene", gene = "PA-X", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, px.e1.end, "0", "+",
        attr_str(gene_id = "PA-X", transcript_id = "unassigned_transcript_2",
                 exception = "ribosomal slippage",
                 gbkey = "CDS", gene = "PA-X", part = "1",
                 product = gene_product("PA-X"), exon_number = "1")),
      gtf_row(seqname, "CDS", px.e2.start, px.e2.end, "0", "+",
        attr_str(gene_id = "PA-X", transcript_id = "unassigned_transcript_2",
                 exception = "ribosomal slippage",
                 gbkey = "CDS", gene = "PA-X", part = "2",
                 product = gene_product("PA-X"), exon_number = "2")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "PA-X", transcript_id = "unassigned_transcript_2",
                 exception = "ribosomal slippage",
                 gbkey = "CDS", gene = "PA-X",
                 product = gene_product("PA-X"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", px.stop.start, px.stop.end, "0", "+",
        attr_str(gene_id = "PA-X", transcript_id = "unassigned_transcript_2",
                 exception = "ribosomal slippage",
                 gbkey = "CDS", gene = "PA-X",
                 product = gene_product("PA-X"), exon_number = "2"))
    )
  }

  # ------------------------------------------------------------------
  # MP: M1 (unspliced) + M2 (spliced mRNA) (canonical len 982)
  # M2 exon1: 1-26, exon2: 715-979, stop: 980-end
  # M1 CDS: 1-756, stop: 757-759
  # ------------------------------------------------------------------
  else if (seg.type == "MP") {

    canonical_len  = 982
    m1.cds.end     = scale_pos(756, canonical_len, seq.len)
    m1.stop.start  = scale_pos(757, canonical_len, seq.len)
    m1.stop.end    = scale_pos(759, canonical_len, seq.len)
    m2.e1.end      = scale_pos(26,  canonical_len, seq.len)
    m2.e2.start    = scale_pos(715, canonical_len, seq.len)
    m2.e2.end      = scale_pos(979, canonical_len, seq.len)
    m2.stop.start  = scale_pos(980, canonical_len, seq.len)
    m2.stop.end    = seq.len

    cat("  Segment:", seqname, "| genes: M1 + M2 (spliced) | len:", seq.len, "\n")
    if (seq.len != canonical_len)
      cat("  Note: scaling MP coordinates from canonical",
          canonical_len, "to", seq.len, "\n")
    check_codon(seq.str, 1,             "start", "M1/M2")
    check_codon(seq.str, m1.stop.start, "stop",  "M1")
    check_codon(seq.str, m2.stop.start, "stop",  "M2")

    rows = c(rows,
      # M2 spliced mRNA (listed first to match reference GTF convention)
      gtf_row(seqname, "gene", 1, seq.len, ".", "+",
        attr_str(gene_id = "M2", transcript_id = "",
                 gbkey = "Gene", gene = "M2", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, m2.e1.end, "0", "+",
        attr_str(gene_id = "M2", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "M2", part = "1",
                 product = gene_product("M2"), exon_number = "1")),
      gtf_row(seqname, "CDS", m2.e2.start, m2.e2.end, "1", "+",
        attr_str(gene_id = "M2", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "M2", part = "2",
                 product = gene_product("M2"), exon_number = "2")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "M2", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "M2",
                 product = gene_product("M2"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", m2.stop.start, m2.stop.end, "0", "+",
        attr_str(gene_id = "M2", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "M2",
                 product = gene_product("M2"), exon_number = "2")),
      # M1 unspliced mRNA
      gtf_row(seqname, "gene", 1, m1.stop.end, ".", "+",
        attr_str(gene_id = "M1", transcript_id = "",
                 gbkey = "Gene", gene = "M1", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, m1.cds.end, "0", "+",
        attr_str(gene_id = "M1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "M1",
                 product = gene_product("M1"), exon_number = "1")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "M1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "M1",
                 product = gene_product("M1"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", m1.stop.start, m1.stop.end, "0", "+",
        attr_str(gene_id = "M1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "M1",
                 product = gene_product("M1"), exon_number = "1"))
    )
  }

  # ------------------------------------------------------------------
  # NS: NS1 (unspliced) + NEP/NS2 (spliced mRNA) (canonical len 838)
  # NEP exon1: 1-30, exon2: 503-835, stop: 836-end
  # NS1 CDS: 1-657, stop: 658-660
  # ------------------------------------------------------------------
  else if (seg.type == "NS") {

    canonical_len  = 838
    ns1.cds.end    = scale_pos(657, canonical_len, seq.len)
    ns1.stop.start = scale_pos(658, canonical_len, seq.len)
    ns1.stop.end   = scale_pos(660, canonical_len, seq.len)
    nep.e1.end     = scale_pos(30,  canonical_len, seq.len)
    nep.e2.start   = scale_pos(503, canonical_len, seq.len)
    nep.e2.end     = scale_pos(835, canonical_len, seq.len)
    nep.stop.start = scale_pos(836, canonical_len, seq.len)
    nep.stop.end   = seq.len

    cat("  Segment:", seqname, "| genes: NS1 + NEP (spliced) | len:", seq.len, "\n")
    if (seq.len != canonical_len)
      cat("  Note: scaling NS coordinates from canonical",
          canonical_len, "to", seq.len, "\n")
    check_codon(seq.str, 1,              "start", "NS1/NEP")
    check_codon(seq.str, ns1.stop.start, "stop",  "NS1")
    check_codon(seq.str, nep.stop.start, "stop",  "NEP")

    rows = c(rows,
      # NEP spliced mRNA (listed first to match reference GTF convention)
      gtf_row(seqname, "gene", 1, seq.len, ".", "+",
        attr_str(gene_id = "NEP", transcript_id = "",
                 gbkey = "Gene", gene = "NEP", gene_biotype = "protein_coding",
                 gene_synonym = "NS2")),
      gtf_row(seqname, "CDS", 1, nep.e1.end, "0", "+",
        attr_str(gene_id = "NEP", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "NEP",
                 note = "nonstructural protein 2", part = "1",
                 product = gene_product("NEP"), exon_number = "1")),
      gtf_row(seqname, "CDS", nep.e2.start, nep.e2.end, "0", "+",
        attr_str(gene_id = "NEP", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "NEP",
                 note = "nonstructural protein 2", part = "2",
                 product = gene_product("NEP"), exon_number = "2")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "NEP", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "NEP",
                 note = "nonstructural protein 2",
                 product = gene_product("NEP"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", nep.stop.start, nep.stop.end, "0", "+",
        attr_str(gene_id = "NEP", transcript_id = "unassigned_transcript_1",
                 gbkey = "CDS", gene = "NEP",
                 note = "nonstructural protein 2",
                 product = gene_product("NEP"), exon_number = "2")),
      # NS1 unspliced mRNA
      gtf_row(seqname, "gene", 1, ns1.stop.end, ".", "+",
        attr_str(gene_id = "NS1", transcript_id = "",
                 gbkey = "Gene", gene = "NS1", gene_biotype = "protein_coding")),
      gtf_row(seqname, "CDS", 1, ns1.cds.end, "0", "+",
        attr_str(gene_id = "NS1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "NS1",
                 product = gene_product("NS1"), exon_number = "1")),
      gtf_row(seqname, "start_codon", 1, 3, "0", "+",
        attr_str(gene_id = "NS1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "NS1",
                 product = gene_product("NS1"), exon_number = "1")),
      gtf_row(seqname, "stop_codon", ns1.stop.start, ns1.stop.end, "0", "+",
        attr_str(gene_id = "NS1", transcript_id = "unassigned_transcript_2",
                 gbkey = "CDS", gene = "NS1",
                 product = gene_product("NS1"), exon_number = "1"))
    )
  }

  return(rows)
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
