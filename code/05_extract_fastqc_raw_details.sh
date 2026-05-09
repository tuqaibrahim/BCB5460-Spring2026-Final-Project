#!/usr/bin/env bash
set -euo pipefail

mkdir -p metadata/fastqc_raw_summaries

for zip in qc/fastqc_raw/*_fastqc.zip; do
    sample=$(basename "$zip" _fastqc.zip)
    out="metadata/fastqc_raw_summaries/${sample}_fastqc_details.txt"

    echo "===== ${sample} =====" > "$out"

    echo -e "\n--- Basic Statistics ---" >> "$out"
    unzip -p "$zip" "${sample}_fastqc/fastqc_data.txt" \
      | awk '/>>Basic Statistics/{flag=1; next} />>END_MODULE/{if(flag){flag=0}} flag' >> "$out"

    echo -e "\n--- Sequence Length Distribution ---" >> "$out"
    unzip -p "$zip" "${sample}_fastqc/fastqc_data.txt" \
      | awk '/>>Sequence Length Distribution/{flag=1; next} />>END_MODULE/{if(flag){flag=0}} flag' >> "$out"

    echo -e "\n--- Adapter Content ---" >> "$out"
    unzip -p "$zip" "${sample}_fastqc/fastqc_data.txt" \
      | awk '/>>Adapter Content/{flag=1; next} />>END_MODULE/{if(flag){flag=0}} flag' >> "$out"

    echo -e "\n--- Overrepresented sequences ---" >> "$out"
    unzip -p "$zip" "${sample}_fastqc/fastqc_data.txt" \
      | awk '/>>Overrepresented sequences/{flag=1; next} />>END_MODULE/{if(flag){flag=0}} flag' >> "$out"

    echo "Wrote $out"
done
