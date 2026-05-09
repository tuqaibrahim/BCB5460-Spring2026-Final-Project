import pandas as pd
import numpy as np
from pathlib import Path

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

pre_cols = ["C1_L-10", "C2_L-10", "C3_L-10"]
post_cols = ["C1_R+3", "C2_R+3", "C3_R+3"]

geo_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/data/geo_processed/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

best_raw_file = Path("counts/from_raw_STAR_featureCounts_gencode_v37_M_O_fraction_noPrimary/lncRNA_counts_s0_gene_name.csv")

geo_de_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/differential_expression/RUVSeq_edgeR_condition_only_LRT_full_results.csv")

raw_de_file = Path("differential_expression_v37_MO_fraction/RUVSeq_edgeR_LRT_candidate_full_results.csv")

outdir = Path("audit_results")
outdir.mkdir(exist_ok=True)

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

# Read count tables
geo = pd.read_csv(geo_file)
raw = pd.read_csv(best_raw_file)

geo = geo.groupby("ID", as_index=False)[sample_cols].sum()
raw = raw.groupby("ID", as_index=False)[sample_cols].sum()

# Read DE tables
geo_de = pd.read_csv(geo_de_file)
raw_de = pd.read_csv(raw_de_file)

geo_de = geo_de[["ID", "logFC", "logCPM", "PValue", "FDR"]].rename(
    columns={
        "logFC": "GEO_DE_logFC",
        "logCPM": "GEO_DE_logCPM",
        "PValue": "GEO_DE_PValue",
        "FDR": "GEO_DE_FDR",
    }
)

raw_de = raw_de[["ID", "logFC", "logCPM", "PValue", "FDR"]].rename(
    columns={
        "logFC": "RAW_DE_logFC",
        "logCPM": "RAW_DE_logCPM",
        "PValue": "RAW_DE_PValue",
        "FDR": "RAW_DE_FDR",
    }
)

# Rename count columns
geo_counts = geo.rename(columns={c: f"{c}_GEO" for c in sample_cols})
raw_counts = raw.rename(columns={c: f"{c}_RAW" for c in sample_cols})

merged = pd.DataFrame({"ID": paper_top10})
merged = merged.merge(geo_counts, on="ID", how="left")
merged = merged.merge(raw_counts, on="ID", how="left")
merged = merged.merge(geo_de, on="ID", how="left")
merged = merged.merge(raw_de, on="ID", how="left")

# Simple paired log2FC helper
def add_simple_fc(df, prefix):
    pre = [f"{c}_{prefix}" for c in pre_cols]
    post = [f"{c}_{prefix}" for c in post_cols]

    df[f"{prefix}_mean_pre"] = df[pre].mean(axis=1)
    df[f"{prefix}_mean_post"] = df[post].mean(axis=1)
    df[f"{prefix}_simple_log2FC_mean_post_over_pre"] = np.log2(
        (df[f"{prefix}_mean_post"] + 0.5) / (df[f"{prefix}_mean_pre"] + 0.5)
    )

    for cpre, cpost, astronaut in zip(pre_cols, post_cols, ["C1", "C2", "C3"]):
        df[f"{prefix}_{astronaut}_paired_log2FC"] = np.log2(
            (df[f"{cpost}_{prefix}"] + 0.5) / (df[f"{cpre}_{prefix}"] + 0.5)
        )

    paired_cols = [f"{prefix}_{a}_paired_log2FC" for a in ["C1", "C2", "C3"]]
    df[f"{prefix}_mean_paired_log2FC"] = df[paired_cols].mean(axis=1)

add_simple_fc(merged, "GEO")
add_simple_fc(merged, "RAW")

# Threshold flags
merged["GEO_paper_like_sig"] = (
    (merged["GEO_DE_logFC"].abs() > 2)
    & (merged["GEO_DE_PValue"] < 0.001)
    & (merged["GEO_DE_FDR"] < 0.05)
)

merged["RAW_paper_like_sig"] = (
    (merged["RAW_DE_logFC"].abs() > 2)
    & (merged["RAW_DE_PValue"] < 0.001)
    & (merged["RAW_DE_FDR"] < 0.05)
)

out = outdir / "paper_top10_GEO_vs_best_raw_audit.csv"
merged.to_csv(out, index=False)

pd.set_option("display.max_columns", 200)
pd.set_option("display.width", 250)

cols_to_show = [
    "ID",
    "GEO_mean_paired_log2FC",
    "RAW_mean_paired_log2FC",
    "GEO_DE_logFC",
    "GEO_DE_PValue",
    "GEO_DE_FDR",
    "GEO_paper_like_sig",
    "RAW_DE_logFC",
    "RAW_DE_PValue",
    "RAW_DE_FDR",
    "RAW_paper_like_sig",
]

print("\nPaper top 10: GEO vs best raw audit")
print(merged[cols_to_show].to_string(index=False))

print("\nSaved:")
print(out)

print("\nCounts for paper top 10:")
count_cols = ["ID"] + [f"{c}_GEO" for c in sample_cols] + [f"{c}_RAW" for c in sample_cols]
print(merged[count_cols].to_string(index=False))
