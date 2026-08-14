#!/usr/bin/env python3
"""Join the depth that gets reported to the depth the callers actually see.

Two different numbers have been called "depth" in this pipeline, and MIN_DEPTH
was stated against the one nothing downstream can see:

  raw            `samtools depth -a -Q 0` -- every aligned base, no quality
                 filter, no overlap handling. This is what depth_profiles/ has
                 always published and what every existing consumer reads.

  (Details of an unpublished run were removed here.)



(Details of an unpublished run were removed here.)





The fix is not to change the floor -- overlap removal is correct, and counting
a fragment's two mates as independent observations is the error it exists to
prevent. The fix is to publish the second number beside the first, so the floor
is stated against a quantity that appears in the output.

This script does the joining. It reads the mpileup on stdin -- the SAME mpileup
command IVAR is given, so the count is the caller's quantity by construction
rather than by a conversion factor -- and writes the raw file back out with the
visible count appended as a fourth column.

(Details of an unpublished run were removed here.)






Note that `-x` is irrelevant to the raw count and must not be used to obtain
it: overlap removal zeroes the mate's base QUALITY without removing the base
from the depth column, so `mpileup -x` and `mpileup` give identical column 4
(verified, 0 disagreements over a full segment). The overlap effect is visible
only after a quality threshold is applied, which is exactly why the two numbers
diverge only once -q enters.

FluLens carries a copy of this counter at tools/reads/depth_band.py, which
measured the band in the first place. Keep the two in step -- they are the same
definition, and if they drift the measurement stops describing the pipeline.

usage:
  samtools mpileup -aa -A -B -d 0 -Q 0 --reference ref.fa in.bam \
    | visible_depth.py --raw sample.depth --min-quality 30 --out sample.4col
"""
import argparse
import sys


def visible(bases, quals, minq):
    """Count bases in a pileup column that survive a -q threshold.

    Walks the base string rather than just counting quality characters,
    because the two are not in step: '^' is followed by a mapping quality and
    '+N'/'-N' by an indel sequence, neither of which consumes a quality, while
    '*' deletion placeholders DO consume one but are not a base observation.
    Counting column 6 alone would include them and overcount beside deletions.
    """
    i = k = keep = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == '^':          # read start, next char is the mapping quality
            i += 2
            continue
        if c == '$':          # read end
            i += 1
            continue
        if c in '+-':         # indel, followed by a length and that many bases
            j = i + 1
            num = ''
            while j < n and bases[j].isdigit():
                num += bases[j]
                j += 1
            i = j + int(num or 0)
            continue
        if c == '*':          # deletion placeholder, consumes a quality
            k += 1
            i += 1
            continue
        if k < len(quals) and (ord(quals[k]) - 33) >= minq:
            keep += 1
        k += 1
        i += 1
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--raw', required=True,
                    help='samtools depth -a output; passed through as columns 1-3')
    ap.add_argument('--min-quality', type=int, default=30,
                    help='base-quality threshold, the same value iVar gets as -q')
    ap.add_argument('--out', required=True, help='output path, or - for stdout')
    args = ap.parse_args()

    # mpileup is read into memory rather than streamed alongside the raw file:
    # both are emitted in BAM-header order and -aa/-a make them the same length,
    # but pairing them line by line would turn any future divergence into a
    # silent off-by-one instead of an error. One flu genome is ~13k positions.
    seen = {}
    for line in sys.stdin:
        f = line.rstrip('\n').split('\t')
        if len(f) < 6:
            continue
        seen[(f[0], f[1])] = visible(f[4], f[5], args.min_quality)

    out = sys.stdout if args.out == '-' else open(args.out, 'w')
    missing = 0
    n = 0
    try:
        with open(args.raw) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line:
                    continue
                f = line.split('\t')
                if len(f) < 3:
                    continue
                vis = seen.get((f[0], f[1]))
                if vis is None:
                    # Not recoverable by guessing: a position the raw pass saw
                    # and the pileup did not means the two are describing
                    # different data, and a fabricated 0 would read as "no
                    # usable coverage here" rather than as the bug it is.
                    missing += 1
                    if missing <= 5:
                        print(f'visible_depth: no pileup row for {f[0]}:{f[1]}',
                              file=sys.stderr)
                    continue
                out.write(f'{line}\t{vis}\n')
                n += 1
    finally:
        if out is not sys.stdout:
            out.close()

    if missing:
        print(f'visible_depth: {missing} position(s) in {args.raw} had no pileup '
              f'row -- the depth pass and the pileup pass disagree about this BAM',
              file=sys.stderr)
        return 1
    print(f'visible_depth: {n} positions, caller-visible depth at -q '
          f'{args.min_quality} written as column 4', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
