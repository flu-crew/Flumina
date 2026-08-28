#!/usr/bin/env python3
"""Which off-subtype HA/NA markers are actually untransferable, and which are
only labelled that way?

FLUMUT_KEEP_MISMATCHED_HA_NA drops every HA/NA marker when the reference's
subtype does not match FluMut's H5/N1 numbering. The case for that is real but
it is an average: FluMut's reference and this run's carry different wild-type
residues at most HA1 positions, so a marker phrased as a deviation from the H5
residue is being read against a protein that never had it.

Per MARKER the answer is already known. flumut_position_map.tsv records, for
every position it places, whether the two references carry the same residue. A
marker sitting on an identical residue is a deviation from the same starting
point, and its H5 provenance is a label rather than a confound.

Reads a finished run and writes one row per dropped HA/NA marker with that
verdict, so the flag stops standing in for a judgement it cannot express.

Rows are counted AFTER reference subtraction — markers_all.tsv minus
reference_markers.tsv. Reading markers_all.tsv alone inflates the answer with
markers this run's own reference carries, which appear in nearly every sample.

usage: hana_transferable.py <run_dir> [out.tsv]
"""
import collections
import os
import re
import sys

MARK = re.compile(r'^([A-Za-z0-9\-]+):([A-Z])(\d+)([A-Z])$')


def rows_of(path, col=1):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > col:
                out.append(f[col])
    return out


def main():
    run = sys.argv[1]
    fm = os.path.join(run, 'variant_analysis', 'flumut')
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(fm, 'hana_transferable.tsv')

    pos = {}
    pmap = os.path.join(fm, 'flumut_position_map.tsv')
    if not os.path.exists(pmap):
        sys.exit(f'no position map at {pmap} -- run FLUMUT_POSITION_MAP first')
    with open(pmap) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 8:
                pos.setdefault(f[0], {})[int(f[2])] = (f[3], f[6], f[7])

    refm = set(rows_of(os.path.join(fm, 'reference_markers.tsv')))
    kept = set(rows_of(os.path.join(fm, 'markers.tsv')))
    counts = collections.Counter(m for m in rows_of(os.path.join(fm, 'markers_all.tsv'))
                                 if m not in refm)

    recs, same = [], 0
    for marker, n in sorted(counts.items(), key=lambda x: -x[1]):
        for tok in (t.strip() for t in marker.split(',')):
            g = MARK.match(tok)
            if not g or g.group(1)[:2] not in ('HA', 'NA'):
                continue
            label, wt, p = g.group(1), g.group(2), int(g.group(3))
            rec = pos.get(label, {}).get(p)
            if rec is None:
                recs.append((marker, tok, label, p, '', '', 'unplaced', n,
                             'yes' if marker in kept else 'no', ''))
                continue
            fm_aa, ref_aa, ident = rec
            if ident == 'yes':
                same += 1
            # The marker's own wild-type should be FluMut's residue. When it is
            # not, the marker is phrased against something other than the
            # reference the map placed, and the verdict is worth less.
            note = '' if wt == fm_aa else f'marker wild-type {wt}'
            recs.append((marker, tok, label, p, fm_aa, ref_aa, ident, n,
                         'yes' if marker in kept else 'no', note))

    with open(out, 'w') as fh:
        fh.write('marker\tcomponent\tlabel\tflumut_pos\tflumut_aa\tref_aa\t'
                 'wt_identical\tsample_rows\tkept\tnote\n')
        for r in recs:
            fh.write('\t'.join(str(x) for x in r) + '\n')

    # One row per COMPONENT, so a compound marker appears twice. Sample-rows are
    # summed over distinct markers instead, or the compound ones count double.
    distinct = {r[0]: r[7] for r in recs}
    print(f'HA/NA after reference subtraction: {len(distinct)} distinct markers, '
          f'{sum(distinct.values())} sample-rows, {len(recs)} components')
    print(f'  components on an IDENTICAL wild-type residue (transferable): {same}')
    print(f'  components on a DIFFERENT wild-type residue                : {len(recs) - same}')
    print(f'written to {out}')


if __name__ == '__main__':
    main()
