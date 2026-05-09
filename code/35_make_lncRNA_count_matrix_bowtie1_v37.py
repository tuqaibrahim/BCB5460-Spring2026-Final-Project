import pandas as pd
import re
from pathlib import Path

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

fc_file = Path("counts/featureCounts_bowtie1_gencode_v37/bowtie1_gencode_v37_featureCounts_s0_M_O_fraction.txt")
gtf_file = Path("reference/gencode_v37/gencode.v37.annotation.gtf")

out_dir = Path("counts/from_bowtie1_featureCounts_gencode_v37_M_O_fraction")
out_dir.mkdir(parents=True, exist_ok=True)

def parse_attr(attr, key):
    m = re.search(rf'{key} "([^"]+)"', attr)
    return m.group(1) if m else None

print("Parsing GENCODE v37 gene annotations...")

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

        records.append({
            "gene_id": parse_attr(attr, "gene_id"),
            "gene_name": parse_attr(attr, "gene_name"),
            "gene_type": parse_attr(attr, "gene_type"),
            "chrom": chrom,
            "start": int(start),
            "end": int(end),
            "strand": strand
        })

anno = pd.DataFrame(records)

print("Gene annotations:", anno.shape)
print("lncRNA genes in annotation:", (anno["gene_type"] == "lncRNA").sum())

print("\nReading Bowtie1 featureCounts table...")
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

merged = fc.merge(
    anno[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]],
    on="gene_id",
    how="left"
)

lnc = merged[merged["gene_type"] == "lncRNA"].copy()

collapsed = (
    lnc.groupby("gene_name", as_index=False)[sample_cols]
    .sum()
    .rename(columns={"gene_name": "ID"})
)

collapsed.insert(1, "Type", "lncRNA")
collapsed["total_count"] = collapsed[sample_cols].sum(axis=1)

count_out = collapsed[collapsed["total_count"] > 0][["ID", "Type"] + sample_cols].copy()

out_lnc = out_dir / "lncRNA_counts_s0_gene_name.csv"
out_info = out_dir / "lncRNA_gene_info_s0.csv"

count_out.to_csv(out_lnc, index=False)
lnc[["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"]].to_csv(out_info, index=False)

print("\nSaved:")
print(out_lnc)
print(out_info)

print("\nShape:", count_out.shape)

print("\nLibrary sizes:")
print(count_out[sample_cols].sum(axis=0))

check_genes = [
    "ZNF436-AS1", "FP671120.7", "AC073140.2", "MIR17HG",
    "MIRLET7BHG", "MIRLET7A1HG", "AD000090.1",
    "LINC00324", "AC024651.2", "LINC02049", "AL732292.2"
]

print("\nSelected genes:")
print(count_out[count_out["ID"].isin(check_genes)].to_string(index=False))
