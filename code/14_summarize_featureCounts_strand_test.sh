#!/usr/bin/env bash
set -euo pipefail

out="metadata/featureCounts_strand_test_summary.tsv"
echo -e "strand\tsample\tassigned\tunassigned_unmapped\tunassigned_multimapping\tunassigned_no_features\tunassigned_ambiguous" > "$out"

for summary in counts/featureCounts_strand_test/*txt.summary; do
    strand=$(basename "$summary" | sed -E 's/.*_s([0-9]+)\.txt\.summary/\1/')

    header=$(head -1 "$summary" | cut -f2-)

    # Loop over sample columns
    ncols=$(awk 'NR==1{print NF}' "$summary")
    for ((i=2; i<=ncols; i++)); do
        sample=$(awk -v i="$i" 'NR==1{print $i}' "$summary")
        sample=$(basename "$sample" _Aligned.sortedByCoord.out.bam)

        assigned=$(awk -v i="$i" '$1=="Assigned"{print $i}' "$summary")
        unmapped=$(awk -v i="$i" '$1=="Unassigned_Unmapped"{print $i}' "$summary")
        multimapping=$(awk -v i="$i" '$1=="Unassigned_MultiMapping"{print $i}' "$summary")
        nofeatures=$(awk -v i="$i" '$1=="Unassigned_NoFeatures"{print $i}' "$summary")
        ambiguous=$(awk -v i="$i" '$1=="Unassigned_Ambiguity"{print $i}' "$summary")

        echo -e "${strand}\t${sample}\t${assigned}\t${unmapped}\t${multimapping}\t${nofeatures}\t${ambiguous}" >> "$out"
    done
done

cat "$out"

