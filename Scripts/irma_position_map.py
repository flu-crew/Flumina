#!/usr/bin/env python3
"""irma_position_map.py

Place IRMA's consensus onto THIS run's reference coordinates, so a downstream
reader can SHOW IRMA's consensus instead of re-deriving one, and can do it by
lookup rather than by assuming a frame.

Why this exists
---------------
IRMA produces a consensus, and it is the consensus FluMut screens. Anything that
paints calls above 50% onto the reference is a SECOND, independent derivation of
the same quantity -- two implementations of one rule, which is the shape of
disagreement this pipeline has been bitten by before. The two are not
interchangeable: where both speak they agree, but IRMA maps to its own
iteratively refined contig and recruits reads BWA soft-clips at the segment
termini, so it carries consensus changes the reference-based callers never see.

Reading IRMA's consensus needs a coordinate frame, and a contig does not
automatically share the reference's. On the 143-sample swine run:

  * 1,084 of 1,123 contigs are exactly reference length,
  * 38 are SHORTER -- and 26 of those 38 by an amount that is not a multiple of
    3, so splicing them through the reference's CDS intervals reads the tail of
    the product out of frame, and
  * one is LONGER: MC-717 A_PA, +1 nt. A single inserted base shifts every
    downstream codon. A truncation at least fails visibly at the end; an
    internal insertion silently returns wrong residues from the insertion
    onward.

So the frame has to be established by alignment, per sample and per segment.
That is a fact to be computed once, here, beside the data -- not re-inferred by
every reader, which is the same conclusion `flumut_position_map.py` reached for
FluMut's numbering, for the same reason.

Depth comes from IRMA's OWN coverage table, not from the BWA pileup
-------------------------------------------------------------------
`min_depth` exists because below it LoFreq and GATK4 report false fixations, and
it is applied to the reference-based alignment. Applying that same floor to an
IRMA consensus call means reading the depth from the alignment that PRODUCED the
call. The two disagree sharply exactly where IRMA is interesting: at `MC-732`
`A_NS` positions 1-5 IRMA has 56-61x where BWA has 2-3x. `<locus>-coverage.txt`
is in the contig's own coordinates and its Consensus column reconstructs the
contig byte-for-byte, so the join is direct.

IRMA has no consensus-depth parameter to point `min_depth` at -- `MIN_TCC` gates
its variants table, not its consensus FASTA -- so the floor is recorded here and
left for the reader to apply. Nothing is masked: 12% of IRMA consensus bases sit
below 100x, and N-masking a FASTA that FluMut screens both manufactures and
suppresses markers.

Output
------
`irma_position_map.tsv` -- one row per sample x product, how the frame was
established and what was found:

    sample product locus status ref_aa_len placed changed ambiguous uncovered thin

    status      identity   contig colinear with the reference, no alignment needed
                aligned    frame established by alignment (truncation or indel)
                absent     no contig for this segment
                no_ref     product's segment is not in the reference FASTA

`irma_consensus_aa.tsv` -- one row per residue that is NOT a plain
reference-matching call at or above the depth floor:

    sample product ref_pos ref_aa irma_aa status irma_depth

    status      change          differs from the reference, depth >= min_depth
                change_thin     differs, but below the floor
                change_nodepth  differs, and depth is UNKNOWN -- no coverage
                                tables were given, so the floor was not applied
                                and this must not be read as passing it
                thin            matches the reference, below the floor
                ambiguous       IRMA has N or a gap here -- NOT a no-change
                uncovered       codon has no placement in the contig at all

A residue absent from that file, in a product whose status is not `absent` or
`no_ref`, matches the reference at or above the floor. `ambiguous` and
`uncovered` are written out rather than left silent on purpose: an absent row
must never be readable as "IRMA agrees", which is the misreading the
zero-coverage marker flag exists to prevent.
"""
import argparse
import os
import re
import sys
from collections import Counter

CODONS = {}
_B, _AA = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for _k, (_i, _j, _l) in enumerate((i, j, l) for i in _B for j in _B for l in _B):
    CODONS[_i + _j + _l] = _AA[_k]


def translate_codon(nt):
    return CODONS.get(nt.upper(), 'X')


def norm_gene(g):
    """Strip a trailing subtype suffix only: HA_H3 -> HA, NA_N2 -> NA.

    Must leave PA-X, PB1-F2, NS1, M1 and M2 untouched -- same rule as makeGTF.R,
    flumut_position_map.py and the viewer, so every consumer agrees on what a
    product is called.
    """
    return re.sub(r'_[HN]?\d+$', '', str(g))


def short_locus(x):
    return re.sub(r'_[A-Z][0-9]+$', '', re.sub(r'^A_', '', str(x or '')))


def read_fasta(path):
    seqs, nm, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                if nm:
                    seqs[nm] = ''.join(buf)
                nm, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if nm:
        seqs[nm] = ''.join(buf)
    return seqs


def read_cds(gtf_dir):
    """CDS intervals per product, exactly as makeGTF.R defines them.

    combined.gtf is skipped and identical intervals de-duplicated, because a
    repeated CDS record would otherwise double a product's length.
    """
    cds, seen = {}, set()
    for fn in sorted(f for f in os.listdir(gtf_dir) if f.endswith('.gtf')):
        if fn.endswith('combined.gtf'):
            continue
        with open(os.path.join(gtf_dir, fn)) as fh:
            for line in fh:
                if not line.strip() or line[0] == '#':
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] != 'CDS':
                    continue
                m = re.search(r'gene_id "([^"]+)"', c[8])
                if not m:
                    continue
                g = norm_gene(m.group(1))
                key = (g, c[0], c[3], c[4])
                if key in seen:
                    continue
                seen.add(key)
                cds.setdefault(g, {'seqname': c[0], 'exons': []})['exons'].append(
                    (int(c[3]), int(c[4])))
    for rec in cds.values():
        rec['exons'].sort()
    return cds


def resolve_locus(name, seqs):
    if name in seqs:
        return name
    for k in seqs:
        if short_locus(k) == short_locus(name):
            return k
    return None


def read_depth(path):
    """IRMA's per-position depth, indexed by CONTIG coordinate (1-based)."""
    depth = []
    try:
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) < 3:
                    continue
                try:
                    depth.append(int(c[2]))
                except ValueError:
                    depth.append(0)
    except OSError:
        return None
    return depth


def make_aligner():
    """Nucleotide aligner, with free end gaps on the REFERENCE.

    A truncated contig has to be free to sit inside the reference without paying
    for the flanks; internal gaps are scored normally, because those are the
    real indels and placing them correctly is the entire point. Mismatch is
    cheap relative to opening a gap, so a divergent stretch is not "explained"
    by a spurious indel.

    `end_insertion_score` is ASSIGNED, never probed -- reading it raises
    ValueError whenever the open and extend end-gap scores differ, as they do
    here, so hasattr() would throw before aligning anything. Same trap as
    flumut_position_map.py.
    """
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = -2.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    aligner.end_insertion_score = 0.0
    return aligner


def colinear_identity(ref_seq, irma_seq):
    same = sum(1 for a, b in zip(ref_seq.upper(), irma_seq.upper()) if a == b)
    return same / max(len(ref_seq), 1)


def map_positions(ref_seq, irma_seq, aligner):
    """ref 0-based index -> contig 0-based index, for positions that align.

    Equal length and high straight identity is taken as colinear without
    aligning: it is 1,084 of 1,123 contigs on a real run, and paying for a
    2,280 x 2,280 alignment to rediscover the identity map is waste. The
    identity threshold is what guards it -- a same-length contig carrying a
    compensating insertion and deletion would fall below it and be aligned
    properly.
    """
    if len(ref_seq) == len(irma_seq) and colinear_identity(ref_seq, irma_seq) >= 0.95:
        return {i: i for i in range(len(ref_seq))}, 'identity'

    best = aligner.align(ref_seq, irma_seq)[0]
    pos = {}
    for (t0, t1), (q0, q1) in zip(*best.aligned):
        for k in range(t1 - t0):
            pos[t0 + k] = q0 + k
    return pos, 'aligned'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference', required=True, help="this run's reference FASTA")
    ap.add_argument('--gtf', required=True, help='reference_gtf directory')
    ap.add_argument('--contigs', required=True,
                    help='IRMA-consensus-contigs directory (<sample>.fasta)')
    ap.add_argument('--irma', default=None,
                    help='IRMA_results directory, for <sample>/tables/<locus>-coverage.txt. '
                         'Without it depths are blank and no residue is called thin')
    ap.add_argument('--min-depth', type=int, default=100,
                    help="depth floor, applied to IRMA's own coverage [100]")
    ap.add_argument('--out', default='irma_position_map.tsv')
    ap.add_argument('--out-aa', default='irma_consensus_aa.tsv')
    args = ap.parse_args()

    if not os.path.isdir(args.contigs):
        sys.stderr.write('irma_position_map: no IRMA consensus contigs; nothing written\n')
        return 0

    seqs = read_fasta(args.reference)
    cds = read_cds(args.gtf)
    if not cds or not seqs:
        sys.stderr.write('irma_position_map: no reference CDS could be read; nothing written\n')
        return 0

    # Reference CDS flattened to a list of 0-based reference indices per product,
    # so a codon is three lookups rather than an interval walk.
    products = {}
    for product, rec in sorted(cds.items()):
        locus = resolve_locus(rec['seqname'], seqs)
        if locus is None:
            products[product] = None
            continue
        idx = []
        for s, e in rec['exons']:
            idx.extend(range(s - 1, e))
        products[product] = {'locus': locus, 'idx': idx}

    try:
        aligner = make_aligner()
    except ImportError:
        sys.stderr.write('irma_position_map: biopython not available; nothing written\n')
        return 0

    map_rows, aa_rows = [], []
    tally = Counter()

    contig_files = sorted(f for f in os.listdir(args.contigs) if f.endswith('.fasta'))
    for fn in contig_files:
        sample = fn[:-6]
        contigs = read_fasta(os.path.join(args.contigs, fn))

        # One alignment per SEGMENT, reused by every product on it -- PA and PA-X
        # share a contig, as do NS1/NEP, M1/M2 and PB1/PB1-F2.
        frames, depths = {}, {}
        for product, rec in products.items():
            if rec is None:
                continue
            locus = rec['locus']
            if locus in frames:
                continue
            con = contigs.get(locus) or contigs.get(resolve_locus(locus, contigs) or '')
            if con is None:
                frames[locus] = (None, 'absent')
                continue
            frames[locus] = map_positions(seqs[locus], con, aligner)
            if args.irma:
                depths[locus] = read_depth(os.path.join(
                    args.irma, sample, 'tables', '%s-coverage.txt' % locus))

        for product, rec in sorted(products.items()):
            if rec is None:
                map_rows.append((sample, product, '', 'no_ref', 0, 0, 0, 0, 0, 0))
                tally['no_ref'] += 1
                continue

            locus = rec['locus']
            pos, status = frames[locus]
            ref_idx = rec['idx']
            n_codons = len(ref_idx) // 3
            ref_seq, con = seqs[locus], contigs.get(locus, '')
            dep = depths.get(locus)

            if pos is None:
                map_rows.append((sample, product, locus, 'absent', n_codons, 0, 0, 0,
                                 n_codons, 0))
                tally['absent'] += 1
                for k in range(n_codons):
                    ref_aa = translate_codon(''.join(ref_seq[i] for i in ref_idx[3 * k:3 * k + 3]))
                    aa_rows.append((sample, product, k + 1, ref_aa, '', 'uncovered', ''))
                continue

            counts = Counter()
            for k in range(n_codons):
                tri = ref_idx[3 * k:3 * k + 3]
                ref_aa = translate_codon(''.join(ref_seq[i] for i in tri))
                mapped = [pos.get(i) for i in tri]

                if any(m is None for m in mapped):
                    counts['uncovered'] += 1
                    aa_rows.append((sample, product, k + 1, ref_aa, '', 'uncovered', ''))
                    continue

                counts['placed'] += 1
                irma_aa = translate_codon(''.join(con[m] for m in mapped))

                d = ''
                thin = False
                if dep:
                    vals = [dep[m] for m in mapped if m < len(dep)]
                    if vals:
                        d = min(vals)
                        thin = d < args.min_depth

                if irma_aa == 'X':
                    counts['ambiguous'] += 1
                    aa_rows.append((sample, product, k + 1, ref_aa, 'X', 'ambiguous', d))
                elif irma_aa != ref_aa:
                    counts['changed'] += 1
                    # An unknown depth is not a passing depth. Without the coverage
                    # tables `change` would read as "at or above the floor" when
                    # nothing was measured — the same absence-of-evidence trap that
                    # `uncovered` exists to keep out of the residue rows.
                    if thin:
                        counts['thin'] += 1
                        state = 'change_thin'
                    else:
                        state = 'change' if d != '' else 'change_nodepth'
                    aa_rows.append((sample, product, k + 1, ref_aa, irma_aa, state, d))
                elif thin:
                    counts['thin'] += 1
                    aa_rows.append((sample, product, k + 1, ref_aa, irma_aa, 'thin', d))

            map_rows.append((sample, product, locus, status, n_codons, counts['placed'],
                             counts['changed'], counts['ambiguous'], counts['uncovered'],
                             counts['thin']))
            tally[status] += 1

    with open(args.out, 'w') as fh:
        fh.write('sample\tproduct\tlocus\tstatus\tref_aa_len\tplaced\tchanged\t'
                 'ambiguous\tuncovered\tthin\n')
        for r in map_rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    with open(args.out_aa, 'w') as fh:
        fh.write('sample\tproduct\tref_pos\tref_aa\tirma_aa\tstatus\tirma_depth\n')
        for r in aa_rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    st = Counter(r[5] for r in aa_rows)
    sys.stderr.write(
        'irma_position_map: %d samples, %d sample x product frames '
        '(%d colinear, %d aligned, %d absent, %d no reference)\n'
        % (len(contig_files), len(map_rows), tally['identity'], tally['aligned'],
           tally['absent'], tally['no_ref']))
    if st['change_nodepth']:
        sys.stderr.write(
            'irma_position_map: NO IRMA coverage tables (--irma not given or unreadable), '
            'so the depth floor was not applied -- %d change_nodepth, %d ambiguous, '
            '%d uncovered\n'
            % (st['change_nodepth'], st['ambiguous'], st['uncovered']))
    else:
        sys.stderr.write(
            'irma_position_map: min_depth %d against IRMA coverage -- '
            '%d change, %d change_thin, %d thin, %d ambiguous, %d uncovered\n'
            % (args.min_depth, st['change'], st['change_thin'], st['thin'],
               st['ambiguous'], st['uncovered']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
