import pandas as pd
from pathlib import Path

infile = Path("counts/from_raw_STAR_featureCounts_gencode_v37_M_O_fraction_noPrimary/lncRNA_counts_s0_gene_name.csv")
outdir = Path("differential_expression_v37_MO_fraction")
outdir.mkdir(exist_ok=True)

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

df = pd.read_csv(infile)

print("Original candidate count matrix shape:", df.shape)
print(df.head())

# Keep only needed columns
df = df[["ID", "Type"] + sample_cols].copy()

# edgeR can technically use non-integer values poorly/incorrectly because counts should be counts.
# featureCounts --fraction gives decimals. GEO table is integer-like.
# For this diagnostic, we will create two versions:
# 1) fractional counts as-is
# 2) rounded counts
#
# The rounded version is closer to what edgeR expects.
frac = df.copy()
rounded = df.copy()
rounded[sample_cols] = rounded[sample_cols].round().astype(int)

# Remove all-zero rows after rounding
rounded["total_count"] = rounded[sample_cols].sum(axis=1)
frac["total_count"] = frac[sample_cols].sum(axis=1)

rounded_nonzero = rounded[rounded["total_count"] > 0].drop(columns=["total_count"]).copy()
frac_nonzero = frac[frac["total_count"] > 0].drop(columns=["total_count"]).copy()

rounded_nonzero.to_csv(outdir / "lncRNA_counts_candidate_rounded.csv", index=False)
frac_nonzero.to_csv(outdir / "lncRNA_counts_candidate_fractional.csv", index=False)

print("\nSaved:")
print(outdir / "lncRNA_counts_candidate_rounded.csv")
print(outdir / "lncRNA_counts_candidate_fractional.csv")

print("\nRounded matrix shape:", rounded_nonzero.shape)
print("Fractional matrix shape:", frac_nonzero.shape)

print("\nRounded library sizes:")
print(rounded_nonzero[sample_cols].sum(axis=0))

print("\nFractional library sizes:")
print(frac_nonzero[sample_cols].sum(axis=0))

paper_top10 = [
    "LINC02049", "AL031767.1", "AC244131.2", "LINC01954", "AL732292.2",
    "AC129926.1", "AC024651.2", "LINC02022", "AC026469.1", "LINC00324"
]

print("\nPaper top 10 genes in rounded candidate:")
print(rounded_nonzero[rounded_nonzero["ID"].isin(paper_top10)].to_string(index=False))
