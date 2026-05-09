import pandas as pd
from pathlib import Path

geo_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/data/geo_processed/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

raw_files = {
    "v37_primary_only": Path("counts/from_raw_STAR_featureCounts_gencode_v37/lncRNA_counts_s0_gene_name.csv"),
    "v37_M": Path("counts/from_raw_STAR_featureCounts_gencode_v37_M/lncRNA_counts_s0_gene_name.csv"),
    "v37_M_fraction": Path("counts/from_raw_STAR_featureCounts_gencode_v37_M_fraction/lncRNA_counts_s0_gene_name.csv"),
}

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

geo = pd.read_csv(geo_file)
geo = geo.groupby("ID", as_index=False)[sample_cols].sum()

out_dir = Path("counts/comparison_to_GEO")
out_dir.mkdir(parents=True, exist_ok=True)

summary_rows = []

paper_top10 = [
    "LINC02049", "AL031767.1", "AC244131.2", "LINC01954", "AL732292.2",
    "AC129926.1", "AC024651.2", "LINC02022", "AC026469.1", "LINC00324"
]

top_geo_genes = [
    "ZNF436-AS1", "FP671120.7", "AC073140.2", "MIR17HG",
    "AD000090.1", "MIRLET7BHG", "MIRLET7A1HG"
]

for mode, raw_file in raw_files.items():
    print(f"\n================ {mode} ================")
    raw = pd.read_csv(raw_file)
    raw = raw.groupby("ID", as_index=False)[sample_cols].sum()

    overlap = sorted(set(raw["ID"]) & set(geo["ID"]))
    print("Raw IDs:", raw["ID"].nunique())
    print("GEO IDs:", geo["ID"].nunique())
    print("Overlap IDs:", len(overlap))
    print("Raw-only:", len(set(raw["ID"]) - set(geo["ID"])))
    print("GEO-only:", len(set(geo["ID"]) - set(raw["ID"])))

    raw2 = raw.rename(columns={c: f"{c}_raw" for c in sample_cols})
    geo2 = geo.rename(columns={c: f"{c}_GEO" for c in sample_cols})
    merged = raw2.merge(geo2, on="ID", how="inner")

    merged["total_raw"] = merged[[f"{c}_raw" for c in sample_cols]].sum(axis=1)
    merged["total_GEO"] = merged[[f"{c}_GEO" for c in sample_cols]].sum(axis=1)
    merged["ratio_raw_over_GEO_total"] = (merged["total_raw"] + 0.5) / (merged["total_GEO"] + 0.5)

    merged.to_csv(out_dir / f"{mode}_raw_vs_GEO_overlap_counts.csv", index=False)

    print("\nLibrary sizes:")
    for c in sample_cols:
        raw_sum = raw[c].sum()
        geo_sum = geo[c].sum()
        print(f"{c}: raw={raw_sum:.2f}, GEO={geo_sum:.2f}, raw/GEO={raw_sum/geo_sum:.3f}")

    print("\nSpearman correlations:")
    spearmans = []
    for c in sample_cols:
        corr = merged[[f"{c}_raw", f"{c}_GEO"]].corr(method="spearman").iloc[0, 1]
        spearmans.append(corr)
        print(f"{c}: {corr:.4f}")

    mean_spearman = sum(spearmans) / len(spearmans)

    summary_rows.append({
        "mode": mode,
        "raw_ids": raw["ID"].nunique(),
        "overlap_ids": len(overlap),
        "mean_spearman": mean_spearman,
        **{f"{c}_raw_sum": raw[c].sum() for c in sample_cols}
    })

    check_genes = paper_top10 + top_geo_genes
    sub = merged[merged["ID"].isin(check_genes)].copy()
    cols = ["ID"] + [f"{c}_raw" for c in sample_cols] + [f"{c}_GEO" for c in sample_cols] + ["total_raw", "total_GEO", "ratio_raw_over_GEO_total"]

    print("\nSelected paper/top-GEO genes:")
    print(sub[cols].sort_values("ID").to_string(index=False))

summary = pd.DataFrame(summary_rows)
summary.to_csv(out_dir / "v37_multimap_mode_comparison_summary.csv", index=False)

print("\n================ SUMMARY ================")
print(summary.to_string(index=False))
