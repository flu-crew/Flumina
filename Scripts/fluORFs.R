#### fluORFs.R
####
#### Single source of truth for influenza A open reading frames.
####
#### Eight segments encode twelve proteins. Four of them are not a simple
#### "translate from nucleotide 1 in frame 1":
####
####   M2      spliced   exon1 1-26,  exon2 715-979  (2 nt carry into exon 2)
####   NEP     spliced   exon1 1-30,  exon2 503-835
####   PA-X    ribosomal slippage, exon1 1-570, exon2 572-697
####   PB1-F2  starts at nucleotide 95, in a different reading frame
####
#### Anything that converts a nucleotide position into an amino-acid position
#### has to go through these intervals. Doing it as ceiling(position/3) is
#### correct ONLY for the eight primary products, and past the end of M1 and NS1
#### it does not fail — it silently returns a plausible number for a codon that
#### does not exist, which is worse.
####
#### This file is sourced by makeGTF.R (which writes these intervals out as GTF),
#### convertVCFtoTable.R and findAAChanges.R so the three cannot drift apart.
#### Base R only, and it derives everything from the reference FASTA — it does
#### NOT read reference_gtf/, which is produced later in the pipeline by the
#### optional SNPGenie step and is therefore not available when the variant
#### table is built.

#############################################
#### Segment identification
#############################################

# Identify the influenza A segment type from a sequence name.
flu_segment_type <- function(name) {
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

# The product a segment is conventionally named for. Variants are reported
# against this one by default, which is what keeps `locus`-keyed joins (the
# curated database, outputSummary.R) behaving exactly as they did before.
flu_primary_product <- function(seg.type) {
  switch(as.character(seg.type),
         "MP" = "M1", "NS" = "NS1",
         as.character(seg.type))
}

# Proportionally scale a canonical coordinate to a different sequence length.
# Same helper makeGTF.R uses; internal ORF coordinates are only meaningful
# relative to the canonical reference they were measured on.
flu_scale_pos <- function(pos, canonical_len, actual_len) {
  round(pos * actual_len / canonical_len)
}

#############################################
#### ORF table
#############################################

# Return the coding intervals for every product on one segment.
#
# Each element of the returned list is:
#   gene     product name (HA, M1, M2, NS1, NEP, PA, PA-X, PB1, PB1-F2, ...)
#   exons    data.frame(start, end) in SEGMENT coordinates, in translation order
#   primary  TRUE for the segment's principal product
#
# Coordinates exclude the stop codon, matching the GTF convention makeGTF.R
# writes. A variant inside a stop codon therefore maps to no product, which is
# reported rather than rounded into the last real codon.
flu_orfs <- function(seq.name, seq.len, seg.type = NULL, seq.str = NULL) {

  if (is.null(seg.type)) seg.type = flu_segment_type(seq.name)
  if (is.na(seg.type))   seg.type = "UNKNOWN"

  ex <- function(starts, ends) data.frame(start = as.integer(starts),
                                          end   = as.integer(ends))

  # ---- single-ORF segments: HA, NA, NP, PB2 (and anything unrecognised) ----
  if (seg.type %in% c("HA", "NA", "NP", "PB2", "UNKNOWN")) {
    gene = if (seg.type == "UNKNOWN") seq.name else seg.type
    return(list(list(gene = gene, exons = ex(1, seq.len - 3), primary = TRUE)))
  }

  # ---- PB1 + PB1-F2 (canonical length 2274) ----
  if (seg.type == "PB1") {
    cl = 2274
    return(list(
      list(gene = "PB1", exons = ex(1, seq.len - 3), primary = TRUE),
      list(gene = "PB1-F2",
           exons = ex(flu_scale_pos(95, cl, seq.len),
                      flu_scale_pos(364, cl, seq.len)),
           primary = FALSE)
    ))
  }

  # ---- PA + PA-X (canonical length 2151) ----
  if (seg.type == "PA") {
    cl = 2151
    return(list(
      list(gene = "PA", exons = ex(1, seq.len - 3), primary = TRUE),
      list(gene = "PA-X",
           exons = ex(c(1,                              flu_scale_pos(572, cl, seq.len)),
                      c(flu_scale_pos(570, cl, seq.len), flu_scale_pos(697, cl, seq.len))),
           primary = FALSE)
    ))
  }

  # ---- MP: M1 unspliced + M2 spliced (canonical length 982) ----
  if (seg.type == "MP") {
    cl = 982
    return(list(
      list(gene = "M1", exons = ex(1, flu_scale_pos(756, cl, seq.len)), primary = TRUE),
      list(gene = "M2",
           exons = ex(c(1,                             flu_scale_pos(715, cl, seq.len)),
                      c(flu_scale_pos(26, cl, seq.len), flu_scale_pos(979, cl, seq.len))),
           primary = FALSE)
    ))
  }

  # ---- NS: NS1 unspliced + NEP/NS2 spliced (canonical length 838) ----
  if (seg.type == "NS") {
    cl = 838
    return(list(
      list(gene = "NS1", exons = ex(1, flu_scale_pos(657, cl, seq.len)), primary = TRUE),
      list(gene = "NEP",
           exons = ex(c(1,                             flu_scale_pos(503, cl, seq.len)),
                      c(flu_scale_pos(30, cl, seq.len), flu_scale_pos(835, cl, seq.len))),
           primary = FALSE)
    ))
  }

  list(list(gene = seq.name, exons = ex(1, seq.len - 3), primary = TRUE))
}

# Wrapper: build the canonical layout, then pull every CDS end back to the first
# in-frame stop when the sequence is available. Callers that have the reference
# should always pass it.
flu_orfs_for <- function(seq.name, seq.str, seg.type = NULL) {
  n = nchar(seq.str)
  orfs = flu_orfs(seq.name, n, seg.type)
  lapply(orfs, function(o) {
    ex = o$exons
    # Run the LAST exon out to the end of the segment before trimming, so the
    # stop search can extend as well as shorten. The canonical end is only a
    # starting guess: NS1 needed extending 657 -> 690 on the H5N1 reference and
    # shortening on others. Trimming alone silently kept the short answer.
    ex$end[nrow(ex)] = n
    o$exons = flu_trim_to_stop(seq.str, ex)
    o
  })
}

#############################################
#### Trimming CDS ends to the real stop codon
#############################################

# The canonical coordinates give the START of each product and its splice
# donor/acceptor sites, which are structurally conserved. They do NOT reliably
# give the END: NS1 and PA-X have strain-variable C-terminal lengths that are
# not proportional to segment length. Both the swine H3N2 and cow H5N1
# references have an 838 nt NS segment, yet NS1 is 219 aa in one and 230 in the
# other; PA-X is likewise 232 aa in one and the full-length 252 in the other.
# Scaling the canonical end truncates the protein and every codon number past
# the cut is then wrong.
#
# So: take the start and the splice sites from the canonical layout, and find
# the end in the sequence itself — the first in-frame stop.
flu_trim_to_stop <- function(seq.str, exons) {
  cds = paste(mapply(function(s, e) substr(seq.str, s, e), exons$start, exons$end),
              collapse = "")
  n = nchar(cds)
  stop.at = NA_integer_
  i = 1L
  while (i + 2L <= n) {
    if (substr(cds, i, i + 2L) %in% c("TAA", "TAG", "TGA")) { stop.at = i; break }
    i = i + 3L
  }
  keep = if (is.na(stop.at)) n - (n %% 3L) else stop.at - 1L
  if (keep <= 0) return(exons)

  out = exons[0, , drop = FALSE]
  used = 0L
  for (k in seq_len(nrow(exons))) {
    len = exons$end[k] - exons$start[k] + 1L
    if (used + len <= keep) {
      out = rbind(out, exons[k, , drop = FALSE]); used = used + len
    } else {
      need = keep - used
      if (need > 0) out = rbind(out, data.frame(start = exons$start[k],
                                                end   = exons$start[k] + need - 1L))
      break
    }
  }
  out
}

#############################################
#### Position mapping
#############################################

# Total coding length of a product, in nucleotides.
flu_cds_length <- function(orf) sum(orf$exons$end - orf$exons$start + 1L)

# Map SEGMENT nucleotide positions onto one product's coding sequence.
# Returns a data.frame with one row per input position:
#   cds_position    1-based index into the spliced CDS, NA if outside it
#   aa_position     codon index, NA if outside
#   codon_position  1/2/3 within the codon, NA if outside
flu_map_positions <- function(orf, pos) {
  pos    = as.integer(pos)
  exons  = orf$exons
  # nucleotides contributed by all preceding exons, so a spliced product
  # numbers straight through its junction
  before = c(0L, cumsum(as.integer(exons$end - exons$start + 1L)))
  cds    = rep(NA_integer_, length(pos))

  for (e in seq_len(nrow(exons))) {
    hit = !is.na(pos) & pos >= exons$start[e] & pos <= exons$end[e]
    if (any(hit)) cds[hit] = before[e] + (pos[hit] - exons$start[e] + 1L)
  }

  data.frame(
    cds_position   = cds,
    aa_position    = ifelse(is.na(cds), NA_integer_, ((cds - 1L) %/% 3L) + 1L),
    codon_position = ifelse(is.na(cds), NA_integer_, ((cds - 1L) %%  3L) + 1L)
  )
}

# The spliced coding sequence of a product, as a plain character string.
# findAAChanges.R translates codons out of this rather than out of the raw
# segment, which is the whole point for M2 / NEP / PA-X.
flu_cds_seq <- function(orf, seq.str) {
  paste(mapply(function(s, e) substr(seq.str, s, e),
               orf$exons$start, orf$exons$end), collapse = "")
}

# Minimal FASTA reader. Kept here so convertVCFtoTable.R kicks off no new
# package dependency just to learn how long each segment is — its only
# requirement stays data.table.
flu_read_fasta <- function(path) {
  ln  = readLines(path)
  ln  = ln[nzchar(ln)]
  idx = grep("^>", ln)
  if (!length(idx)) stop("No FASTA records found in: ", path)
  nm  = sub("\\s.*$", "", sub("^>", "", ln[idx]))
  st  = idx + 1L
  en  = c(idx[-1] - 1L, length(ln))
  sq  = mapply(function(a, b) if (a > b) "" else paste(ln[a:b], collapse = ""), st, en)
  stats::setNames(toupper(sq), nm)
}

#############################################
#### Annotating a variant table
#############################################

# Attach product / amino-acid annotation to a table of variant calls.
#
# `dat` needs a segment-name column and a nucleotide-position column. The
# returned table gains:
#   product          the gene the position codes for
#   product_primary  TRUE when that is the segment's principal product
#   cds_position     1-based index into that product's spliced CDS
#   aa_position      codon index within that product
#   codon_position   1/2/3 within the codon
#
# By default a position that codes in two products yields TWO rows, one per
# product — that is the only honest representation of overlapping reading
# frames, and it is what makes M2, NEP, PA-X and PB1-F2 visible at all. With
# primary.only = TRUE the table keeps exactly one row per input row, which is
# what callers that key on (locus, position) need.
#
# Positions in no coding region (UTR, stop codon, or a spliced-out intron) are
# dropped when expanding, and carry NA under primary.only. They are never given
# a made-up codon number.
flu_annotate_positions <- function(dat, ref, locus.col = "locus",
                                   pos.col = "position", primary.only = FALSE) {

  if (!nrow(dat)) return(cbind(dat, product = character(), product_primary = logical(),
                               cds_position = integer(), aa_position = integer(),
                               codon_position = integer()))

  loci = unique(dat[[locus.col]])
  miss = loci[!loci %in% names(ref)]
  if (length(miss))
    warning("Segment(s) absent from the reference, left unannotated: ",
            paste(miss, collapse = ", "))

  out = list()
  for (lc in loci) {
    sub = dat[dat[[locus.col]] == lc, , drop = FALSE]
    if (!lc %in% names(ref)) {
      sub$product = NA_character_; sub$product_primary = NA
      sub$cds_position = NA_integer_; sub$aa_position = NA_integer_
      sub$codon_position = NA_integer_
      out[[length(out) + 1L]] = sub
      next
    }
    orfs = flu_orfs_for(lc, ref[[lc]])
    if (primary.only) orfs = Filter(function(o) o$primary, orfs)

    for (o in orfs) {
      m    = flu_map_positions(o, sub[[pos.col]])
      keep = if (primary.only) rep(TRUE, nrow(sub)) else !is.na(m$aa_position)
      if (!any(keep)) next
      piece = sub[keep, , drop = FALSE]
      piece$product         = o$gene
      piece$product_primary = o$primary
      piece$cds_position    = m$cds_position[keep]
      piece$aa_position     = m$aa_position[keep]
      piece$codon_position  = m$codon_position[keep]
      out[[length(out) + 1L]] = piece
    }
  }
  res = do.call(rbind, out)

  # Put the new columns beside what they describe, and do it HERE rather than in
  # the caller so a table annotated on the way through convertVCFtoTable.R and
  # one back-filled later by findAAChanges.R come out with identical layout.
  lead = c("method", "sample", "locus", "product", "product_primary",
           "position", "reference", "alternative", "quality", "depth",
           "map_quality", "allele_frequency",
           "cds_position", "aa_position", "codon_position")
  res[, c(lead[lead %in% colnames(res)], setdiff(colnames(res), lead))]
}

# Every product of every segment in a reference, as one flat data.frame.
# Convenience for callers that want to look up by (locus, product).
flu_orf_table <- function(seq.names, seq.lens) {
  out = list()
  for (i in seq_along(seq.names)) {
    for (o in flu_orfs(seq.names[i], seq.lens[i])) {
      out[[length(out) + 1L]] = data.frame(
        locus      = seq.names[i],
        product    = o$gene,
        primary    = o$primary,
        n_exons    = nrow(o$exons),
        cds_length = flu_cds_length(o),
        aa_length  = flu_cds_length(o) %/% 3L,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}
