import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# INPUT / OUTPUT
# ==========================
counts_file = "/lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv"

out_png = "/lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline/lncRNA_10panel_CPM_plot.png"
out_pdf = "/lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline/lncRNA_10panel_CPM_plot.pdf"

# ==========================
# CHOOSE 5 UP + 5 DOWN GENES
# ==========================
# Example set based on the paper-style panel you showed
up_genes = [
"LINC02049",
"AL031767.1",
"AC244131.2",
"LINC01954",
"AL732292.2"
]

down_genes = [
"AC129926.1",
"AC024651.2",
"LINC02022",
"AC026469.1",
"LINC00324"
]

genes_to_plot = up_genes + down_genes

# ==========================
# READ COUNTS FILE
# ==========================
df = pd.read_csv(counts_file)

print("Columns in file:")
print(df.columns.tolist())

# Rename ID column to gene_id to make code cleaner
df = df.rename(columns={"ID": "gene_id"})

if "gene_id" not in df.columns:
raise ValueError("Could not find gene column. Expected column named 'ID'.")

# ==========================
# DEFINE SAMPLE COLUMNS
# ==========================
L10_samples = [c for c in df.columns if c.endswith("_L-10")]
R3_samples = [c for c in df.columns if c.endswith("_R+3")]
count_cols = L10_samples + R3_samples

print("L-10 samples:", L10_samples)
print("R+3 samples:", R3_samples)

if len(L10_samples) == 0 or len(R3_samples) == 0:
raise ValueError("Could not detect sample columns for L-10 and R+3.")

# Make sure counts are numeric
for col in count_cols:
df[col] = pd.to_numeric(df[col], errors="coerce")

# ==========================
# CALCULATE CPM USING ALL GENES
# ==========================
library_sizes = df[count_cols].sum(axis=0)
print("\nLibrary sizes:")
print(library_sizes)

cpm_all = df[count_cols].div(library_sizes, axis=1) * 1e6

# Build CPM dataframe
cpm_df = pd.concat([df[["gene_id", "Type"]], cpm_all], axis=1)

# ==========================
# SUBSET TO TARGET GENES
# ==========================
plot_df = cpm_df[cpm_df["gene_id"].isin(genes_to_plot)].copy()

missing = [g for g in genes_to_plot if g not in plot_df["gene_id"].values]
if missing:
print("\nMissing genes from counts file:")
print(missing)

# Keep requested order
plot_df["gene_id"] = pd.Categorical(plot_df["gene_id"], categories=genes_to_plot, ordered=True)
plot_df = plot_df.sort_values("gene_id")

# ==========================
# CONVERT TO LONG FORMAT
# ==========================
long_rows = []

for _, row in plot_df.iterrows():
gene = row["gene_id"]

for sample in L10_samples:
long_rows.append({
"gene_id": gene,
"group": "L-10",
"sample": sample,
"CPM": row[sample]
})

for sample in R3_samples:
long_rows.append({
"gene_id": gene,
"group": "R+3",
"sample": sample,
"CPM": row[sample]
})

long_df = pd.DataFrame(long_rows)

# ==========================
# SUMMARY: mean + SEM
# ==========================
summary = (
long_df
.groupby(["gene_id", "group"])
.agg(
mean_CPM=("CPM", "mean"),
sd_CPM=("CPM", "std"),
n=("CPM", "count")
)
.reset_index()
)

summary["sem_CPM"] = summary["sd_CPM"] / np.sqrt(summary["n"])

# ==========================
# PLOT
# ==========================
fig, axes = plt.subplots(2, 5, figsize=(16, 7))
axes = axes.flatten()

bar_colors = ["#b7d7ea", "#4a86c5"] # light blue, blue

for i, gene in enumerate(genes_to_plot):
ax = axes[i]

sub = summary[summary["gene_id"] == gene].copy()

if sub.empty:
ax.axis("off")
ax.set_title(f"{gene}\n(not found)", fontsize=11, color="gray")
continue

# order groups
sub["group"] = pd.Categorical(sub["group"], categories=["L-10", "R+3"], ordered=True)
sub = sub.sort_values("group")

means = sub["mean_CPM"].values
errors = sub["sem_CPM"].values

ax.bar(
["L-10", "R+3"],
means,
yerr=errors,
capsize=5,
color=bar_colors,
edgecolor="none"
)

title_color = "darkgreen" if gene in up_genes else "firebrick"
ax.set_title(gene, color=title_color, fontsize=13, fontweight="bold")
ax.set_ylabel("CPM", fontsize=11)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=10)

# remove extra spines for cleaner look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Optional panel label like paper
fig.text(0.02, 0.96, "E", fontsize=24, fontweight="bold")

plt.tight_layout(rect=[0.03, 0.02, 1, 0.95])
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("\nSaved figure to:")
print(out_png)
print(out_pdf)