#!/usr/bin/env python3
"""flumut_position_map.py

Emit the mapping from FluMut's marker numbering onto THIS run's reference
coordinates, so downstream readers join markers to variants by lookup instead of
inferring an offset.

Why this exists
---------------
A FluMut marker is written `HA1-5:G224S` -- protein HA1, residue 224, in the
numbering of FluMut's own reference. Nothing in the marker string says where
residue 224 falls in the reference a given run was mapped against. Downstream
this was being recovered by searching for a constant shift that put every
marker's residue onto the translated reference, which fails in two ways:

  * it needs several markers per protein to pin a shift down, so small proteins
    (M2, PB1-F2, NS-2) never resolve at all; and
  * a constant shift is only correct when the two proteins are colinear. HA and
    NA have diverged by INDELS between subtypes, so no single offset is right
    across the protein -- on an H3N2 run the best shift for HA1 explained 8 of
    14 markers, which is neither a fit nor a clean failure.

Both problems disappear here. FluMut ships its reference sequences and their CDS
annotations in flumut_db.sqlite, so the mapping is a fact to be read and aligned,
not a parameter to be estimated. Aligning protein-to-protein handles the indels,
which is what makes HA and NA work across subtypes rather than only on H5N1.

Output: flumut_position_map.tsv, one row per FluMut residue that aligns.

    label       marker prefix as it appears in results (`HA1-5`, `NS-1`, `PB2`)
    protein     FluMut protein the label belongs to, per the database
    flumut_pos  residue number as written in the marker
    flumut_aa   residue FluMut's reference carries there
    product     product name in THIS run's GTF
    ref_pos     codon number in THIS run's CDS  <- the join key
    ref_aa      residue THIS run's reference carries there
    identical   whether the two references agree at that residue

A position absent from the file did not align and must not be placed. That is
deliberate: an absent position is recoverable, a confidently wrong one is not.
"""
import argparse
import os
import re
import sqlite3
import sys

CODONS = {}
_B, _AA = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for _k, (_i, _j, _l) in enumerate((i, j, l) for i in _B for j in _B for l in _B):
    CODONS[_i + _j + _l] = _AA[_k]


def translate(cds):
    return ''.join(CODONS.get(cds[i:i + 3].upper(), 'X') for i in range(0, len(cds) - 2, 3))


def norm_gene(g):
    """Strip a trailing subtype suffix only: HA_H3 -> HA, NA_N2 -> NA.

    Must leave PA-X, PB1-F2, NS1, M1 and M2 untouched -- same rule as makeGTF.R
    and the viewer, so all three agree on what a product is called.
    """
    return re.sub(r'_[HN]?\d+$', '', str(g))


def short_locus(x):
    return re.sub(r'_[A-Z][0-9]+$', '', re.sub(r'^A_', '', str(x or '')))


# FluMut protein -> product name in our GTF. FluMut splits HA into its two mature
# subunits and numbers each from 1; our GTF has one HA ORF, so both map onto it
# and the alignment is what places them.
PRODUCT = {'HA1': 'HA', 'HA2': 'HA', 'NA': 'NA', 'NS-1': 'NS1', 'NS-2': 'NEP',
           'M1': 'M1', 'M2': 'M2', 'NP': 'NP', 'PA': 'PA', 'PA-X': 'PA-X',
           'PB1': 'PB1', 'PB1-F2': 'PB1-F2', 'PB2': 'PB2'}


def find_db(explicit):
    if explicit:
        return explicit
    try:
        import flumutdb
        p = os.path.join(os.path.dirname(flumutdb.__file__), 'flumut_db.sqlite')
        if os.path.exists(p):
            return p
    except ImportError:
        pass
    return None


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


def read_our_proteins(reference_fa, gtf_dir):
    """Translate this run's CDS exactly as makeGTF.R defines them.

    combined.gtf is skipped and identical intervals are de-duplicated, because a
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

    seqs = read_fasta(reference_fa)
    prot = {}
    for g, rec in cds.items():
        rec['exons'].sort()
        seq = seqs.get(rec['seqname'])
        if seq is None:
            for k, v in seqs.items():
                if short_locus(k) == short_locus(rec['seqname']):
                    seq = v
                    break
        if seq is None:
            continue
        prot[g] = translate(''.join(seq[s - 1:e] for s, e in rec['exons']))
    return prot


def read_flumut_proteins(db):
    """Splice and translate each protein from FluMut's own reference sequences."""
    con = sqlite3.connect(db)
    refs = {name: seq for name, seq in con.execute('select name, sequence from "references"')}
    ann = {}
    for start, end, protein, refname in con.execute(
            'select start, end, protein_name, reference_name from annotations'):
        ann.setdefault(protein, {'ref': refname, 'exons': []})['exons'].append((start, end))

    labels = {}
    for name, protein in con.execute('select name, protein_name from mutations'):
        labels.setdefault(name.split(':')[0], protein)

    version = None
    try:
        version = '.'.join(str(x) for x in list(
            con.execute('select major, minor from db_version'))[0])
    except (sqlite3.Error, IndexError):
        pass
    con.close()

    prot = {}
    for protein, rec in ann.items():
        seq = refs.get(rec['ref'])
        if seq is None:
            continue
        rec['exons'].sort()
        prot[protein] = translate(''.join(seq[s - 1:e] for s, e in rec['exons']))
    return prot, labels, version


def align(our_seq, fm_seq):
    """Align FluMut's protein onto ours, returning [(our_idx0, fm_idx0), ...].

    End gaps on OUR sequence are free: HA1 and HA2 are subunits of a protein our
    GTF holds as one ORF, so the query has to be free to sit inside the target
    without paying for the flanks. Internal gaps are scored normally -- those are
    the real indels between subtypes, and they are the whole reason this is an
    alignment rather than an offset.
    """
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = 'global'
    # Free end gaps on OUR sequence (Biopython >= 1.86 spelling; the pinned image
    # supplies it). Assign, never probe: READING this attribute raises ValueError
    # whenever the open and extend end-gap scores differ, as they do here, so
    # hasattr() would throw before aligning anything.
    aligner.end_insertion_score = 0.0

    best = aligner.align(our_seq, fm_seq)[0]
    pairs = []
    for (t0, t1), (q0, q1) in zip(*best.aligned):
        for k in range(t1 - t0):
            pairs.append((t0 + k, q0 + k))
    return pairs, best.score


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference', required=True, help='this run\'s reference FASTA')
    ap.add_argument('--gtf', required=True, help='reference_gtf directory')
    ap.add_argument('--db', default=None, help='flumut_db.sqlite (default: from flumutdb)')
    ap.add_argument('--out', default='flumut_position_map.tsv')
    args = ap.parse_args()

    db = find_db(args.db)
    if not db:
        sys.stderr.write('flumut_position_map: flumut database not found; '
                         'no position map written\n')
        return 0

    our = read_our_proteins(args.reference, args.gtf)
    fm, labels, dbver = read_flumut_proteins(db)
    if not our:
        sys.stderr.write('flumut_position_map: no reference proteins could be translated; '
                         'no position map written\n')
        return 0

    rows, summary = [], []
    for label, protein in sorted(labels.items()):
        product = PRODUCT.get(protein)
        fm_seq = fm.get(protein)
        our_seq = our.get(product) if product else None
        if not fm_seq or not our_seq:
            summary.append(f'  {label:<8} {protein:<8} -> '
                           f'{"no FluMut protein" if not fm_seq else f"no {product} in this reference"}')
            continue

        pairs, _ = align(our_seq, fm_seq)
        ident = sum(1 for t, q in pairs if our_seq[t] == fm_seq[q])
        for t, q in pairs:
            rows.append((label, protein, q + 1, fm_seq[q], product, t + 1, our_seq[t],
                         'yes' if our_seq[t] == fm_seq[q] else 'no'))

        offsets = {t - q for t, q in pairs}
        shape = (f'offset {list(offsets)[0]:+d}' if len(offsets) == 1
                 else f'{len(offsets)} blocks (indels)')
        summary.append(f'  {label:<8} {protein:<8} -> {product:<7} '
                       f'{len(pairs):>4}/{len(fm_seq)} aligned, '
                       f'{100 * ident / max(len(pairs), 1):5.1f}% identical, {shape}')

    with open(args.out, 'w') as fh:
        fh.write('label\tprotein\tflumut_pos\tflumut_aa\tproduct\tref_pos\tref_aa\tidentical\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    sys.stderr.write(f'flumut_position_map: flumut db {dbver or "unknown"}, '
                     f'{len(rows)} positions across {len(set(r[0] for r in rows))} labels\n')
    sys.stderr.write('\n'.join(summary) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
