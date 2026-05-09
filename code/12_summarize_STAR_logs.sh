#!/usr/bin/env bash
set -euo pipefail

out="metadata/STAR_alignment_summary.tsv"
echo -e "sample\tinput_reads\tunique_mapped_percent\tmulti_mapped_reads\ttoo_short_percent" > "$out"

for log in star_alignments/*Log.final.out; do
    sample=$(basename "$log" _Log.final.out)

    input=$(grep "Number of input reads" "$log" | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')
    unique=$(grep "Uniquely mapped reads %" "$log" | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')
    multi=$(grep "Number of reads mapped to multiple loci" "$log" | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')
    short=$(grep "% of reads unmapped: too short" "$log" | awk -F'|' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}')

    echo -e "${sample}\t${input}\t${unique}\t${multi}\t${short}" >> "$out"
done

cat "$out"
