import pandas as pd
import re
from pathlib import Path

# --------------------------
# Input files
# --------------------------
fc_file = Path("counts/featureCounts_strand_test/gencode_v38_featureCounts_s0.txt")
gtf_file = Path("reference/gencode_v38/gencode.v38.primary_assembly.annotation.gtf")

out_dir = Path("counts/from_raw_STAR_featureCounts")
out_dir.mkdir(parents=True, exist_ok=True)

out_all = out_dir / "all_gene_counts_s0.csv"
out_lnc = out_dir / "lncRNA_counts_s0_gene_name.csv"
out_gene_info = out_dir / "lncRNA_gene_info_s0.csv"

# --------------------------
# Parse GENCODE GTF attributes
# --------------------------
def parse_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else None

print("Parsing GTF gene annotations...")

records = []
with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attr = parts

        if feature != "gene":
            continue

        gene_id = parse_attr(attr, "gene_id")
        gene_name = parse_attr(attr, "gene_name")
        gene_type = parse_attr(attr, "gene_type")

        records.append({
            "gene_id": gene_id,
            "gene_id_no_version": gene_id.split(".")[0] if gene_id else None,
            "gene_name": gene_name,
            "gene_type": gene_type,
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand
        })

anno = pd.DataFrame(records)

print("Number of gene annotations:", len(anno))
print("\nGene type counts:")
print(anno["gene_type"].value_counts().head(30))

# --------------------------
# Read featureCounts table
# --------------------------
print("\nReading featureCounts table...")

fc = pd.read_csv(fc_file, sep="\t", comment="#")

print("featureCounts shape:", fc.shape)
print("featureCounts columns:")
print(fc.columns.tolist())

# Rename sample columns from BAM path names to clean sample names
rename = {}
for col in fc.columns:
    if col.endswith("_Aligned.sortedByCoord.out.bam"):
        sample = Path(col).name.replace("_Aligned.sortedByCoord.out.bam", "")
        rename[col] = sample

fc = fc.rename(columns=rename)

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

missing = [c for c in sample_cols if c not in fc.columns]
if missing:
    raise ValueError(f"Missing expected sample columns: {missing}")

# featureCounts Geneid should be GENCODE gene_id
fc = fc.rename(columns={"Geneid": "gene_id"})

# Save all gene counts
all_counts = fc[["gene_id"] + sample_cols].copy()
all_counts.to_csv(out_all, index=False)

# --------------------------
# Join annotation
# --------------------------
merged = fc.merge(
    anno[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]],
    on="gene_id",
    how="left"
)

print("\nAfter joining annotation:")
print(merged["gene_type"].value_counts(dropna=False).head(30))

# --------------------------
# Keep lncRNA genes
# --------------------------
# GENCODE uses gene_type == "lncRNA" for the broad lncRNA class.
lnc = merged[merged["gene_type"] == "lncRNA"].copy()

print("\nNumber of lncRNA genes:", len(lnc))

# Check missing names
print("Missing gene_name among lncRNAs:", lnc["gene_name"].isna().sum())

# --------------------------
# Use gene_name as ID because GEO table used names like LINC02049, AC024651.2
# Some gene_names can duplicate, so collapse by summing counts.
# --------------------------
before = len(lnc)
dup_names = lnc["gene_name"].duplicated().sum()

print("lncRNA rows before collapsing duplicate gene_names:", before)
print("Number of duplicated lncRNA gene_names:", dup_names)

collapsed = (
    lnc.groupby("gene_name", as_index=False)[sample_cols]
    .sum()
    .rename(columns={"gene_name": "ID"})
)

collapsed.insert(1, "Type", "lncRNA")

# Remove all-zero rows
collapsed["total_count"] = collapsed[sample_cols].sum(axis=1)
n_zero = (collapsed["total_count"] == 0).sum()
print("All-zero lncRNA rows:", n_zero)

collapsed_nonzero = collapsed[collapsed["total_count"] > 0].copy()

# Save count matrix
count_out = collapsed_nonzero[["ID", "Type"] + sample_cols].copy()
count_out.to_csv(out_lnc, index=False)

# Save gene info
gene_info = lnc[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]].copy()
gene_info.to_csv(out_gene_info, index=False)

print("\nSaved:")
print(out_all)
print(out_lnc)
print(out_gene_info)

print("\nlncRNA count matrix shape:", count_out.shape)
print(count_out.head())

print("\nLibrary sizes from raw-derived lncRNA matrix:")
print(count_out[sample_cols].sum(axis=0))
