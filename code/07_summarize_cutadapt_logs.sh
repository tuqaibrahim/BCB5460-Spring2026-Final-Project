#!/usr/bin/env bash
set -euo pipefail

out="metadata/cutadapt_summary.tsv"
echo -e "sample\ttotal_reads\treads_written\tpercent_written\treads_too_short\tadapter1_illumina\tadapter2_polyA\tadapter3_polyG" > "$out"

for log in cutadapt_logs/*_cutadapt.log; do
    sample=$(basename "$log" _cutadapt.log)

    total=$(grep "Total reads processed:" "$log" | awk '{print $4}' | tr -d ',')
    written_line=$(grep "Reads written (passing filters):" "$log")
    written=$(echo "$written_line" | awk '{print $5}' | tr -d ',')
    percent=$(echo "$written_line" | grep -oP '\(\K[0-9.]+(?=%\))' || true)

    too_short=$(grep "Reads that were too short:" "$log" | awk '{print $6}' | tr -d ',' || echo "NA")

    # Count adapter sections if available
    adapter1=$(grep -A 2 "Sequence: AGATCGGAAGAGC" "$log" | grep "Trimmed:" | awk '{print $2}' | tr -d ',' || echo "NA")
    adapter2=$(grep -A 2 "Sequence: AAAAAAAAAA" "$log" | grep "Trimmed:" | awk '{print $2}' | tr -d ',' || echo "NA")
    adapter3=$(grep -A 2 "Sequence: GGGGGGGGGG" "$log" | grep "Trimmed:" | awk '{print $2}' | tr -d ',' || echo "NA")

    echo -e "${sample}\t${total}\t${written}\t${percent}\t${too_short}\t${adapter1}\t${adapter2}\t${adapter3}" >> "$out"
done

cat "$out"

