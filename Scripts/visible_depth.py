#!/usr/bin/env python3
"""Join the depth that gets reported to the depth the callers actually see.

Two different numbers get called "depth" here:

  raw            `samtools depth -a -Q 0` -- every aligned base. What
                 depth_profiles/ publishes and every existing consumer reads.
  caller-visible `samtools mpileup` with overlapping mate bases zeroed, then
                 -q applied. What iVar's -m floor is actually tested against,
                 and well below the raw count.

MIN_DEPTH was stated against the second while only the first was published, so
positions could clear the floor on paper and be evaluated by no caller. This
script reads the mpileup on stdin -- the SAME command IVAR is given, so the
count is the caller's quantity rather than a conversion -- and appends it to
the raw file as a fourth column.

Columns 1-3 pass through byte-for-byte. mpileup's own column 4 is NOT used for
them: it counts deletion placeholders where `samtools depth` does not, and
swapping one for the other would move an already-published number.

`-x` will not give you the raw count either -- overlap removal zeroes the
mate's base QUALITY without removing it from the depth column, so the two
diverge only once -q is applied.

FluLens tools/reads/depth_band.py carries a copy of this counter. Keep them in
step; if they drift, that measurement stops describing this pipeline.

Full measurement and reasoning in the Flumina HANDOFF.md.

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

    # Read into memory rather than streamed alongside the raw file: pairing the
    # two line by line would turn any future divergence into a silent
    # off-by-one instead of an error. One flu genome is ~13k positions.
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
                    # Never fabricate a 0 here: it would read as "no usable
                    # coverage" rather than as the bug it is.
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
