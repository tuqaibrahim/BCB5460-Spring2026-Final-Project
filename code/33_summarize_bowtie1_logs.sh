#!/usr/bin/env bash
set -euo pipefail

out="metadata/bowtie1_alignment_summary.tsv"
echo -e "sample\tlog_file" > "$out"

for log in bowtie1_alignments/*_bowtie1.log; do
    sample=$(basename "$log" _bowtie1.log)
    echo -e "${sample}\t${log}" >> "$out"
done

cat "$out"

echo
echo "===== Bowtie1 log contents ====="
for log in bowtie1_alignments/*_bowtie1.log; do
    echo
    echo "----- $log -----"
    cat "$log"
done
