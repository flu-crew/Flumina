#!/usr/bin/env bash
# Compare Snakemake-era output against Nextflow output on the same input.
#
# A Nextflow port that runs is not a port that is correct. This diffs the actual
# variant calls, ignoring the things that legitimately differ between runs
# (file paths, timestamps, command lines, tool invocation order in headers).
#
#   ./dev/compare_pipelines.sh  old_output_dir  new_output_dir
#
# Exit 0 = variant calls identical. Exit 1 = they differ (details printed).

set -uo pipefail

OLD="${1:?usage: compare_pipelines.sh OLD_DIR NEW_DIR}"
NEW="${2:?usage: compare_pipelines.sh OLD_DIR NEW_DIR}"

for d in "$OLD" "$NEW"; do
    [ -d "$d" ] || { echo "not a directory: $d" >&2; exit 2; }
done

# Strip VCF header lines and any column that encodes a path or a run timestamp.
# What remains is the biology: CHROM POS REF ALT QUAL FILTER and the INFO field.
normalise_vcf() {
    grep -v '^##' "$1" 2>/dev/null \
      | awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1,$2,$4,$5,$6,$7,$8}' \
      | sort
}

status=0
total=0; same=0; differ=0; missing=0

echo "comparing variant calls"
echo "  old: $OLD"
echo "  new: $NEW"
echo

# Iterate the union of samples present in either tree
samples=$( { ls "$OLD/vcf_files" 2>/dev/null; ls "$NEW/vcf_files" 2>/dev/null; } | sort -u )
[ -z "$samples" ] && { echo "no vcf_files/ in either directory" >&2; exit 2; }

for sample in $samples; do
    for vcf in gatk4-filtered-snps.vcf lofreq-called-variants.vcf gatk4-unfiltered-genotypes.vcf; do
        a="$OLD/vcf_files/$sample/$vcf"
        b="$NEW/vcf_files/$sample/$vcf"
        total=$((total+1))

        if [ ! -f "$a" ] || [ ! -f "$b" ]; then
            printf '  %-8s %-34s %s\n' "MISSING" "$sample" "$vcf"
            [ -f "$a" ] || echo "             absent in old"
            [ -f "$b" ] || echo "             absent in new"
            missing=$((missing+1)); status=1
            continue
        fi

        if diff -q <(normalise_vcf "$a") <(normalise_vcf "$b") >/dev/null 2>&1; then
            na=$(normalise_vcf "$a" | wc -l | tr -d ' ')
            printf '  %-8s %-34s %-34s %s calls\n' "OK" "$sample" "$vcf" "$na"
            same=$((same+1))
        else
            printf '  %-8s %-34s %s\n' "DIFFER" "$sample" "$vcf"
            # show a bounded sample of the difference rather than dumping everything
            diff <(normalise_vcf "$a") <(normalise_vcf "$b") \
              | head -12 | sed 's/^/             /'
            only_a=$(comm -23 <(normalise_vcf "$a") <(normalise_vcf "$b") | wc -l | tr -d ' ')
            only_b=$(comm -13 <(normalise_vcf "$a") <(normalise_vcf "$b") | wc -l | tr -d ' ')
            echo "             $only_a call(s) only in old, $only_b only in new"
            differ=$((differ+1)); status=1
        fi
    done
done

echo
echo "----------------------------------------------------------"
printf 'compared %s files: %s identical, %s differing, %s missing\n' \
       "$total" "$same" "$differ" "$missing"
if [ "$status" -eq 0 ]; then
    echo "RESULT: variant calls are identical — port is faithful"
else
    echo "RESULT: outputs differ — investigate before trusting the port"
fi
exit "$status"
