import pandas as pd
from pathlib import Path

outdir = Path("audit_results")
outdir.mkdir(exist_ok=True)

rows = [
    {
        "pipeline": "GEO processed matrix",
        "annotation": "GEO-provided processed lncRNA table",
        "aligner": "unknown from processed table",
        "featureCounts_mode": "unknown",
        "overlap_with_GEO": "reference",
        "mean_spearman_vs_GEO": "reference",
        "DE_result": "GEO matrix contains paper-like signal for some top genes",
        "main_conclusion": "Use as reference processed matrix."
    },
    {
        "pipeline": "STAR + GENCODE v38 primary + featureCounts primary-only",
        "annotation": "GENCODE v38 primary",
        "aligner": "STAR",
        "featureCounts_mode": "-s 0 --primary",
        "overlap_with_GEO": 3421,
        "mean_spearman_vs_GEO": "~0.44",
        "DE_result": "Does not reproduce paper genes",
        "main_conclusion": "Poor ID overlap; many paper/GEO IDs absent from v38 annotation."
    },
    {
        "pipeline": "STAR + GENCODE v37 + featureCounts primary-only",
        "annotation": "GENCODE v37",
        "aligner": "STAR",
        "featureCounts_mode": "-s 0 --primary",
        "overlap_with_GEO": 9225,
        "mean_spearman_vs_GEO": 0.382,
        "DE_result": "Does not reproduce paper genes",
        "main_conclusion": "Annotation overlap improves, but counts still mismatch."
    },
    {
        "pipeline": "STAR + GENCODE v37 + featureCounts -M --primary",
        "annotation": "GENCODE v37",
        "aligner": "STAR",
        "featureCounts_mode": "-s 0 -M --primary",
        "overlap_with_GEO": 10300,
        "mean_spearman_vs_GEO": 0.442,
        "DE_result": "Not selected for final DE",
        "main_conclusion": "Counting multimappers improves overlap but still misses key GEO-high genes."
    },
    {
        "pipeline": "STAR + GENCODE v37 + featureCounts -M -O --fraction",
        "annotation": "GENCODE v37",
        "aligner": "STAR",
        "featureCounts_mode": "-s 0 -M -O --fraction, no --primary",
        "overlap_with_GEO": 11038,
        "mean_spearman_vs_GEO": 0.504,
        "DE_result": "34 significant genes; none of paper top 10 significant",
        "main_conclusion": "Best raw-derived reconstruction, but still does not reproduce paper DE genes."
    },
    {
        "pipeline": "Bowtie1 + GENCODE v37 + featureCounts -M -O --fraction",
        "annotation": "GENCODE v37",
        "aligner": "Bowtie1",
        "featureCounts_mode": "-s 0 -M -O --fraction",
        "overlap_with_GEO": 10595,
        "mean_spearman_vs_GEO": 0.431,
        "DE_result": "Not run for DE because count match was worse than STAR",
        "main_conclusion": "Bowtie1 did not improve global match to GEO."
    },
]

df = pd.DataFrame(rows)

out_csv = outdir / "pipeline_audit_summary.csv"
out_tsv = outdir / "pipeline_audit_summary.tsv"

df.to_csv(out_csv, index=False)
df.to_csv(out_tsv, sep="\t", index=False)

print(df.to_string(index=False))
print("\nSaved:")
print(out_csv)
print(out_tsv)
