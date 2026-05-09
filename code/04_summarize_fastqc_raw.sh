#!/usr/bin/env bash
set -euo pipefail

out="metadata/fastqc_raw_summary_status.tsv"
echo -e "sample\tstatus\tmodule\tfilename" > "$out"

for zip in qc/fastqc_raw/*_fastqc.zip; do
    sample=$(basename "$zip" _fastqc.zip)

    unzip -p "$zip" "${sample}_fastqc/summary.txt" \
      | awk -v s="$sample" 'BEGIN{OFS="\t"} {status=$1; filename=$NF; module=$0; sub("^[^\t]*\t","",module); sub("\t[^\t]*$","",module); print s,status,module,filename}' \
      >> "$out"
done

cat "$out"
