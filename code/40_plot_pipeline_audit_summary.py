import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

infile = Path("audit_results/pipeline_audit_summary.csv")
outdir = Path("audit_results")
outdir.mkdir(exist_ok=True)

df = pd.read_csv(infile)

# Keep only numeric pipelines, excluding GEO reference row
plot_df = df[df["overlap_with_GEO"] != "reference"].copy()
plot_df["overlap_with_GEO"] = pd.to_numeric(plot_df["overlap_with_GEO"])
plot_df["mean_spearman_vs_GEO"] = (
    plot_df["mean_spearman_vs_GEO"]
    .astype(str)
    .str.replace("~", "", regex=False)
    .astype(float)
)

# Short labels for plotting
label_map = {
    "STAR + GENCODE v38 primary + featureCounts primary-only": "STAR v38\nprimary-only",
    "STAR + GENCODE v37 + featureCounts primary-only": "STAR v37\nprimary-only",
    "STAR + GENCODE v37 + featureCounts -M --primary": "STAR v37\n-M primary",
    "STAR + GENCODE v37 + featureCounts -M -O --fraction": "STAR v37\n-M -O fraction",
    "Bowtie1 + GENCODE v37 + featureCounts -M -O --fraction": "Bowtie1 v37\n-M -O fraction",
}
plot_df["label"] = plot_df["pipeline"].map(label_map)

# Sort in logical order
plot_df = plot_df.set_index("label").loc[
    [
        "STAR v38\nprimary-only",
        "STAR v37\nprimary-only",
        "STAR v37\n-M primary",
        "STAR v37\n-M -O fraction",
        "Bowtie1 v37\n-M -O fraction",
    ]
].reset_index()

# Figure 1: overlap with GEO
plt.figure(figsize=(9, 5))
plt.bar(plot_df["label"], plot_df["overlap_with_GEO"])
plt.ylabel("Number of lncRNA IDs overlapping GEO")
plt.xlabel("Pipeline")
plt.title("Pipeline audit: lncRNA ID overlap with GEO processed matrix")
plt.xticks(rotation=35, ha="right")
plt.tight_layout()
plt.savefig(outdir / "pipeline_audit_overlap_with_GEO.png", dpi=300)
plt.savefig(outdir / "pipeline_audit_overlap_with_GEO.pdf")
plt.close()

# Figure 2: Spearman correlation
plt.figure(figsize=(9, 5))
plt.bar(plot_df["label"], plot_df["mean_spearman_vs_GEO"])
plt.ylabel("Mean Spearman correlation vs GEO")
plt.xlabel("Pipeline")
plt.title("Pipeline audit: count similarity to GEO processed matrix")
plt.xticks(rotation=35, ha="right")
plt.ylim(0, 0.6)
plt.tight_layout()
plt.savefig(outdir / "pipeline_audit_spearman_vs_GEO.png", dpi=300)
plt.savefig(outdir / "pipeline_audit_spearman_vs_GEO.pdf")
plt.close()

print("Saved:")
print(outdir / "pipeline_audit_overlap_with_GEO.png")
print(outdir / "pipeline_audit_overlap_with_GEO.pdf")
print(outdir / "pipeline_audit_spearman_vs_GEO.png")
print(outdir / "pipeline_audit_spearman_vs_GEO.pdf")

print("\nPlot data:")
print(plot_df[["label", "overlap_with_GEO", "mean_spearman_vs_GEO"]].to_string(index=False))
