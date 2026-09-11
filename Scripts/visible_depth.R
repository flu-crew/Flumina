#!/usr/bin/env Rscript

#### visible_depth.R
####
#### Join reported depth with the depth visible to the callers.
####
#### Two different numbers are called "depth" here:
####   raw            `samtools depth -a -Q 0` -- every aligned base. What
####                  depth_profiles/ publishes and every existing consumer reads.
####   caller-visible `samtools mpileup` with overlapping mate bases zeroed, then
####                  -q applied. What iVar's -m floor is tested against, and well
####                  below the raw count.
####
#### MIN_DEPTH is defined against the caller-visible value while only the raw
#### value was published, so a position could clear the floor on paper and be
#### evaluated by no caller. This script reads the mpileup on stdin -- the same
#### command iVar is given -- and appends the caller-visible count to the raw file
#### as a fourth column.
####
#### Columns 1-3 pass through unchanged. mpileup's own column 4 is not used for
#### them: it counts deletion placeholders where `samtools depth` does not, and
#### swapping one for the other would move an already-published number.
####
#### FluLens tools/reads/depth_band.py carries a copy of this counter. Keep them
#### in step.
####
#### usage:
####   samtools mpileup -aa -A -B -d 0 -Q 0 --reference ref.fa in.bam \
####     | Rscript visible_depth.R --raw sample.depth --min-quality 30 --out sample.4col

# Remove the '+N'/'-N' indel markers from a pileup base string, along with the N
# bases each one introduces. The inserted bases are nucleotides, so they never
# form another '[+-][0-9]+' match, and the ranges do not overlap.
remove.indels = function(s) {
  loc = gregexpr("[+-][0-9]+", s, perl = TRUE)[[1]]
  if (loc[1] == -1L) return(s)
  marker.len = attr(loc, "match.length")
  markers = regmatches(s, gregexpr("[+-][0-9]+", s, perl = TRUE))[[1]]
  counts = as.integer(sub("^[+-]", "", markers))
  ns = nchar(s)
  del.start = loc
  del.end = pmin(loc + marker.len - 1L + counts, ns)

  # Keep the substrings between the deleted ranges.
  kept = character(0)
  pos = 1L
  for (k in seq_along(del.start)) {
    if (del.start[k] > pos) kept = c(kept, substr(s, pos, del.start[k] - 1L))
    pos = max(pos, del.end[k] + 1L)
  }
  if (pos <= ns) kept = c(kept, substr(s, pos, ns))
  paste0(kept, collapse = "")
}

# Count the bases in one pileup column whose quality is at or above the -q floor.
#
# The base string and the quality string are not in step: '^' is followed by a
# mapping-quality character and '+N'/'-N' by an indel sequence, neither of which
# consumes a quality, while '*' deletion placeholders consume one but are not a
# base observation. Remove the read-start, read-end, and indel markers, which
# leaves one character per quality character (a base or a '*'). Then count the
# bases, but not the '*', whose quality passes the floor.
visible = function(bases, quals, minq) {
  b = gsub("\\^.", "", bases, perl = TRUE)   # read start plus its mapping quality
  b = gsub("$", "", b, fixed = TRUE)         # read end
  b = remove.indels(b)

  b.raw = charToRaw(b)
  q.scores = as.integer(charToRaw(quals)) - 33L
  n = min(length(b.raw), length(q.scores))
  if (n == 0L) return(0L)
  sum(b.raw[seq_len(n)] != charToRaw("*") & q.scores[seq_len(n)] >= minq)
}

# Read a `--flag value` pair from the argument vector.
args = commandArgs(trailingOnly = TRUE)
get.flag = function(name, default = NULL) {
  hit = which(args == name)
  if (length(hit) == 0L) return(default)
  i = hit[length(hit)]
  if (i == length(args)) {
    cat("Error: ", name, " needs a value\n", sep = "", file = stderr())
    quit(status = 1)
  }
  args[i + 1L]
}

raw.path = get.flag("--raw")
out.path = get.flag("--out")
min.quality = suppressWarnings(as.integer(get.flag("--min-quality", "30")))

if (is.null(raw.path) || is.null(out.path)) {
  cat("Usage: samtools mpileup ... | Rscript visible_depth.R --raw <depth> ",
      "--min-quality <int> --out <path|->\n", sep = "", file = stderr())
  quit(status = 1)
}
if (is.na(min.quality)) {
  cat("Error: --min-quality must be an integer\n", file = stderr())
  quit(status = 1)
}

# Read the pileup from stdin into a lookup keyed by chromosome and position. A
# flu genome is about 13,000 positions, so holding it in memory is cheap.
seen = new.env(parent = emptyenv())
con = file("stdin", "r")
pileup = readLines(con)
close(con)
for (line in pileup) {
  f = strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(f) < 6L) next
  assign(paste(f[1], f[2], sep = "\t"), visible(f[5], f[6], min.quality), envir = seen)
}

# Pass the raw depth file through, appending the caller-visible count as column 4.
out.con = if (out.path == "-") stdout() else file(out.path, "w")
raw = readLines(raw.path)
missing = 0L
written = 0L
for (line in raw) {
  if (!nzchar(line)) next
  f = strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(f) < 3L) next
  key = paste(f[1], f[2], sep = "\t")
  if (!exists(key, envir = seen, inherits = FALSE)) {
    # Do not fabricate zero: it would read as no usable coverage rather than as
    # a missing pileup row.
    missing = missing + 1L
    if (missing <= 5L) {
      cat("visible_depth: no pileup row for ", f[1], ":", f[2], "\n",
          sep = "", file = stderr())
    }
    next
  }
  cat(line, "\t", get(key, envir = seen, inherits = FALSE), "\n",
      sep = "", file = out.con)
  written = written + 1L
}
if (out.path != "-") close(out.con)

if (missing > 0L) {
  cat("visible_depth: ", missing, " position(s) in ", raw.path,
      " had no pileup row -- the depth pass and the pileup pass disagree about",
      " this BAM\n", sep = "", file = stderr())
  quit(status = 1)
}
cat("visible_depth: ", written, " positions, caller-visible depth at -q ",
    min.quality, " written as column 4\n", sep = "", file = stderr())
quit(status = 0)
