#!/usr/bin/env bash
set -euo pipefail

out="metadata/merged_fastq_read_counts.tsv"
echo -e "sample\tfastq\treads" > "$out"

for fq in raw_merged_fastq/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    lines=$(zcat "$fq" | wc -l)
    reads=$((lines / 4))
    echo -e "${sample}\t${fq}\t${reads}" >> "$out"
done

cat "$out"
