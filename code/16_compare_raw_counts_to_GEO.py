import pandas as pd
from pathlib import Path

# --------------------------
# Files
# --------------------------
raw_file = Path("counts/from_raw_STAR_featureCounts/lncRNA_counts_s0_gene_name.csv")

# Adjust this path if needed
geo_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/data/geo_processed/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

out_dir = Path("counts/comparison_to_GEO")
out_dir.mkdir(parents=True, exist_ok=True)

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

# --------------------------
# Read tables
# --------------------------
raw = pd.read_csv(raw_file)
geo = pd.read_csv(geo_file)

print("Raw-derived shape:", raw.shape)
print("GEO shape:", geo.shape)

# Keep lncRNA rows only if Type exists
if "Type" in geo.columns:
    print("\nGEO Type counts:")
    print(geo["Type"].value_counts())

# Collapse GEO duplicate IDs if any
geo_collapsed = geo.groupby("ID", as_index=False)[sample_cols].sum()
raw_collapsed = raw.groupby("ID", as_index=False)[sample_cols].sum()

print("\nRaw-derived unique IDs:", raw_collapsed["ID"].nunique())
print("GEO unique IDs:", geo_collapsed["ID"].nunique())

overlap = sorted(set(raw_collapsed["ID"]) & set(geo_collapsed["ID"]))
print("Overlap IDs:", len(overlap))

print("Raw-only IDs:", len(set(raw_collapsed["ID"]) - set(geo_collapsed["ID"])))
print("GEO-only IDs:", len(set(geo_collapsed["ID"]) - set(raw_collapsed["ID"])))

# --------------------------
# Library sizes
# --------------------------
print("\nLibrary sizes:")
lib = pd.DataFrame({
    "sample": sample_cols,
    "raw_derived_lncRNA_sum": [raw_collapsed[c].sum() for c in sample_cols],
    "GEO_lncRNA_sum": [geo_collapsed[c].sum() for c in sample_cols],
})
lib["raw_over_GEO"] = lib["raw_derived_lncRNA_sum"] / lib["GEO_lncRNA_sum"]
print(lib.to_string(index=False))
lib.to_csv(out_dir / "library_size_comparison.csv", index=False)

# --------------------------
# Merge overlap
# --------------------------
raw2 = raw_collapsed.rename(columns={c: f"{c}_raw" for c in sample_cols})
geo2 = geo_collapsed.rename(columns={c: f"{c}_GEO" for c in sample_cols})

merged = raw2.merge(geo2, on="ID", how="inner")

# Total counts across samples
merged["total_raw"] = merged[[f"{c}_raw" for c in sample_cols]].sum(axis=1)
merged["total_GEO"] = merged[[f"{c}_GEO" for c in sample_cols]].sum(axis=1)
merged["ratio_raw_over_GEO_total"] = (merged["total_raw"] + 0.5) / (merged["total_GEO"] + 0.5)

# Save full overlap
merged.to_csv(out_dir / "raw_vs_GEO_overlap_counts.csv", index=False)

# --------------------------
# Correlations per sample
# --------------------------
print("\nSpearman correlations raw vs GEO by sample:")
corr_rows = []
for c in sample_cols:
    corr = merged[[f"{c}_raw", f"{c}_GEO"]].corr(method="spearman").iloc[0, 1]
    corr_rows.append({"sample": c, "spearman": corr})
    print(c, corr)

pd.DataFrame(corr_rows).to_csv(out_dir / "raw_vs_GEO_spearman_by_sample.csv", index=False)

# --------------------------
# Paper genes
# --------------------------
paper_top10 = [
    "LINC02049",
    "AL031767.1",
    "AC244131.2",
    "LINC01954",
    "AL732292.2",
    "AC129926.1",
    "AC024651.2",
    "LINC02022",
    "AC026469.1",
    "LINC00324"
]

paper = merged[merged["ID"].isin(paper_top10)].copy()
print("\nPaper top 10 raw vs GEO counts:")
cols = ["ID"] + [f"{c}_raw" for c in sample_cols] + [f"{c}_GEO" for c in sample_cols] + ["total_raw", "total_GEO", "ratio_raw_over_GEO_total"]
print(paper[cols].sort_values("ID").to_string(index=False))
paper[cols].to_csv(out_dir / "paper_top10_raw_vs_GEO_counts.csv", index=False)

missing_paper = sorted(set(paper_top10) - set(merged["ID"]))
print("\nPaper top 10 missing from overlap:")
print(missing_paper)

# --------------------------
# Top raw and GEO genes
# --------------------------
print("\nTop 20 GEO lncRNAs by total count:")
print(merged.sort_values("total_GEO", ascending=False)[["ID", "total_raw", "total_GEO", "ratio_raw_over_GEO_total"]].head(20).to_string(index=False))

print("\nTop 20 raw-derived lncRNAs by total count:")
print(merged.sort_values("total_raw", ascending=False)[["ID", "total_raw", "total_GEO", "ratio_raw_over_GEO_total"]].head(20).to_string(index=False))
