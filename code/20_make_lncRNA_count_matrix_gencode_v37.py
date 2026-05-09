import pandas as pd
import re
from pathlib import Path

fc_file = Path("counts/featureCounts_gencode_v37/gencode_v37_featureCounts_s0.txt")
gtf_file = Path("reference/gencode_v37/gencode.v37.annotation.gtf")

out_dir = Path("counts/from_raw_STAR_featureCounts_gencode_v37")
out_dir.mkdir(parents=True, exist_ok=True)

out_all = out_dir / "all_gene_counts_s0.csv"
out_lnc = out_dir / "lncRNA_counts_s0_gene_name.csv"
out_gene_info = out_dir / "lncRNA_gene_info_s0.csv"

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

def parse_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else None

print("Parsing GENCODE v37 GTF gene annotations...")

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

print("\nReading featureCounts table...")
fc = pd.read_csv(fc_file, sep="\t", comment="#")

rename = {}
for col in fc.columns:
    if col.endswith("_Aligned.sortedByCoord.out.bam"):
        sample = Path(col).name.replace("_Aligned.sortedByCoord.out.bam", "")
        rename[col] = sample

fc = fc.rename(columns=rename)
fc = fc.rename(columns={"Geneid": "gene_id"})

missing = [c for c in sample_cols if c not in fc.columns]
if missing:
    raise ValueError(f"Missing expected sample columns: {missing}")

all_counts = fc[["gene_id"] + sample_cols].copy()
all_counts.to_csv(out_all, index=False)

merged = fc.merge(
    anno[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]],
    on="gene_id",
    how="left"
)

print("\nAfter joining annotation:")
print(merged["gene_type"].value_counts(dropna=False).head(30))

lnc = merged[merged["gene_type"] == "lncRNA"].copy()

print("\nNumber of lncRNA genes:", len(lnc))
print("Missing gene_name among lncRNAs:", lnc["gene_name"].isna().sum())

print("lncRNA rows before collapsing duplicate gene_names:", len(lnc))
print("Number of duplicated lncRNA gene_names:", lnc["gene_name"].duplicated().sum())

collapsed = (
    lnc.groupby("gene_name", as_index=False)[sample_cols]
    .sum()
    .rename(columns={"gene_name": "ID"})
)

collapsed.insert(1, "Type", "lncRNA")
collapsed["total_count"] = collapsed[sample_cols].sum(axis=1)

print("All-zero lncRNA rows:", (collapsed["total_count"] == 0).sum())

count_out = collapsed[collapsed["total_count"] > 0][["ID", "Type"] + sample_cols].copy()
count_out.to_csv(out_lnc, index=False)

gene_info = lnc[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]].copy()
gene_info.to_csv(out_gene_info, index=False)

print("\nSaved:")
print(out_all)
print(out_lnc)
print(out_gene_info)

print("\nlncRNA count matrix shape:", count_out.shape)
print(count_out.head())

print("\nLibrary sizes from raw-derived v37 lncRNA matrix:")
print(count_out[sample_cols].sum(axis=0))

paper_genes = [
    "AC024651.2", "AC244131.2", "AL031767.1", "AL732292.2",
    "AC129926.1", "AC026469.1", "LINC00324", "LINC02049"
]

print("\nPaper genes in v37 raw-derived count matrix:")
print(count_out[count_out["ID"].isin(paper_genes)].to_string(index=False))
