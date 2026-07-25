#!/usr/bin/env python3
"""
Compare two Flumina output directories (e.g. Snakemake vs Nextflow).

This goes further than diffing files. The specific risk when porting is that
per-sample results get attached to the WRONG sample — every file present, every
count plausible, but the labels shuffled. A naive diff of sorted variant tables
would not catch that. So this checks set equality *per sample*, and separately
tests whether one pipeline's sample-A calls turn up under sample B in the other.

    ./dev/compare_outputs.py OLD_DIR NEW_DIR

Exit 0 only if per-sample call sets are identical in both directions.
"""
import csv
import sys
from pathlib import Path


def read_vcf_calls(path):
    """(CHROM, POS, REF, ALT) for each record. Headers and INFO ignored."""
    calls = set()
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) >= 5:
                    calls.add((f[0], f[1], f[3], f[4]))
    except FileNotFoundError:
        return None
    return calls


def read_variant_table(path):
    """{sample: {(locus, position, reference, alternative, method)}}"""
    by_sample = {}
    try:
        with open(path, newline='') as fh:
            for row in csv.DictReader(fh):
                key = (row.get('locus'), row.get('position'),
                       row.get('reference'), row.get('alternative'),
                       row.get('method'))
                by_sample.setdefault(row.get('sample'), set()).add(key)
    except FileNotFoundError:
        return None
    return by_sample


def banner(t):
    print(f"\n{t}\n" + "-" * len(t))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    old, new = Path(sys.argv[1]), Path(sys.argv[2])
    for d in (old, new):
        if not d.is_dir():
            sys.exit(f"not a directory: {d}")

    failures = []

    # ---------------------------------------------------------------- VCFs
    banner("per-sample VCF calls")
    samples = sorted({p.name for d in (old, new)
                      for p in (d / 'vcf_files').glob('*') if p.is_dir()})
    if not samples:
        sys.exit("no vcf_files/<sample>/ directories found in either tree")

    for vcf_name in ('gatk4-filtered-snps.vcf', 'lofreq-called-variants.vcf'):
        print(f"\n  {vcf_name}")
        for s in samples:
            a = read_vcf_calls(old / 'vcf_files' / s / vcf_name)
            b = read_vcf_calls(new / 'vcf_files' / s / vcf_name)
            if a is None or b is None:
                where = 'old' if a is None else 'new'
                print(f"    MISSING  {s[:44]:<44} (absent in {where})")
                failures.append(f"{s}/{vcf_name} missing in {where}")
                continue
            if a == b:
                print(f"    OK       {s[:44]:<44} {len(a)} calls")
            else:
                print(f"    DIFFER   {s[:44]:<44} old={len(a)} new={len(b)} "
                      f"only_old={len(a - b)} only_new={len(b - a)}")
                for c in sorted(a - b)[:3]:
                    print(f"               only in old: {c}")
                for c in sorted(b - a)[:3]:
                    print(f"               only in new: {c}")
                failures.append(f"{s}/{vcf_name} differs")

    # ------------------------------------------------- variant table by sample
    banner("variant-table.csv, per sample")
    ta = read_variant_table(old / 'variant_analysis' / 'variant-table.csv')
    tb = read_variant_table(new / 'variant_analysis' / 'variant-table.csv')

    if ta is None or tb is None:
        print("  variant-table.csv missing on one side — skipping table checks")
        failures.append("variant-table.csv missing")
    else:
        if set(ta) != set(tb):
            print(f"  SAMPLE SETS DIFFER\n    only old: {sorted(set(ta)-set(tb))}"
                  f"\n    only new: {sorted(set(tb)-set(ta))}")
            failures.append("sample sets differ")
        for s in sorted(set(ta) & set(tb)):
            if ta[s] == tb[s]:
                print(f"  OK       {s[:44]:<44} {len(ta[s])} calls")
            else:
                print(f"  DIFFER   {s[:44]:<44} old={len(ta[s])} new={len(tb[s])}")
                failures.append(f"table rows differ for {s}")

        # ------------------------------------------------ shuffle detection
        # If labels were swapped, sample X's calls in `old` would match some
        # OTHER sample's calls in `new`. Check every cross pair.
        banner("shuffle check (cross-sample matching)")
        shuffled = False
        for sx in sorted(set(ta) & set(tb)):
            if ta[sx] == tb[sx]:
                continue
            for sy in sorted(tb):
                if sy != sx and ta[sx] and ta[sx] == tb[sy]:
                    print(f"  SHUFFLED: old '{sx}' == new '{sy}'")
                    failures.append(f"labels shuffled: {sx} -> {sy}")
                    shuffled = True
        if not shuffled:
            print("  no cross-sample matches — labels are not shuffled")

    # ------------------------------------------------------------- verdict
    banner("verdict")
    if failures:
        print(f"  {len(failures)} problem(s):")
        for f in failures:
            print(f"    - {f}")
        return 1
    print("  identical: per-sample call sets match in both pipelines")
    return 0


if __name__ == '__main__':
    sys.exit(main())
