import pandas as pd
from pathlib import Path

geo_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/data/geo_processed/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

raw_files = {
    "STAR_v37_M_O_fraction_noPrimary": Path("counts/from_raw_STAR_featureCounts_gencode_v37_M_O_fraction_noPrimary/lncRNA_counts_s0_gene_name.csv"),
    "Bowtie1_v37_M_O_fraction": Path("counts/from_bowtie1_featureCounts_gencode_v37_M_O_fraction/lncRNA_counts_s0_gene_name.csv"),
}

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

geo = pd.read_csv(geo_file)
geo = geo.groupby("ID", as_index=False)[sample_cols].sum()

selected_genes = [
    "ZNF436-AS1", "FP671120.7", "AC073140.2", "MIR17HG",
    "AD000090.1", "MIRLET7BHG", "MIRLET7A1HG",
    "LINC00324", "AC024651.2", "LINC02049", "AL732292.2"
]

summary_rows = []

for mode, raw_file in raw_files.items():
    print(f"\n================ {mode} ================")
    raw = pd.read_csv(raw_file)
    raw = raw.groupby("ID", as_index=False)[sample_cols].sum()

    raw2 = raw.rename(columns={c: f"{c}_raw" for c in sample_cols})
    geo2 = geo.rename(columns={c: f"{c}_GEO" for c in sample_cols})

    merged = raw2.merge(geo2, on="ID", how="inner")

    merged["total_raw"] = merged[[f"{c}_raw" for c in sample_cols]].sum(axis=1)
    merged["total_GEO"] = merged[[f"{c}_GEO" for c in sample_cols]].sum(axis=1)
    merged["ratio_raw_over_GEO_total"] = (merged["total_raw"] + 0.5) / (merged["total_GEO"] + 0.5)

    spearmans = []
    for c in sample_cols:
        corr = merged[[f"{c}_raw", f"{c}_GEO"]].corr(method="spearman").iloc[0, 1]
        spearmans.append(corr)

    mean_spearman = sum(spearmans) / len(spearmans)

    print("Raw IDs:", raw["ID"].nunique())
    print("GEO IDs:", geo["ID"].nunique())
    print("Overlap IDs:", len(merged))
    print("Mean Spearman:", round(mean_spearman, 4))

    print("\nLibrary raw/GEO ratios:")
    for c in sample_cols:
        print(f"{c}: {raw[c].sum() / geo[c].sum():.3f}")

    sub = merged[merged["ID"].isin(selected_genes)].copy()
    cols = (
        ["ID"] +
        [f"{c}_raw" for c in sample_cols] +
        [f"{c}_GEO" for c in sample_cols] +
        ["total_raw", "total_GEO", "ratio_raw_over_GEO_total"]
    )

    print("\nSelected genes:")
    print(sub[cols].sort_values("ID").to_string(index=False))

    summary_rows.append({
        "mode": mode,
        "raw_ids": raw["ID"].nunique(),
        "overlap_ids": len(merged),
        "mean_spearman": mean_spearman,
        **{f"{c}_raw_sum": raw[c].sum() for c in sample_cols},
        **{f"{c}_raw_over_GEO": raw[c].sum() / geo[c].sum() for c in sample_cols},
    })

summary = pd.DataFrame(summary_rows).sort_values("mean_spearman", ascending=False)

print("\n================ SUMMARY ================")
print(summary.to_string(index=False))

Path("counts/comparison_to_GEO").mkdir(parents=True, exist_ok=True)
summary.to_csv("counts/comparison_to_GEO/bowtie1_vs_star_best_vs_GEO_summary.csv", index=False)
