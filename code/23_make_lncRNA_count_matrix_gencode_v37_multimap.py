import pandas as pd
import re
from pathlib import Path

sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]

gtf_file = Path("reference/gencode_v37/gencode.v37.annotation.gtf")

modes = {
    "M": Path("counts/featureCounts_gencode_v37_multimap/gencode_v37_featureCounts_s0_M.txt"),
    "M_fraction": Path("counts/featureCounts_gencode_v37_multimap/gencode_v37_featureCounts_s0_M_fraction.txt"),
}

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

for mode, fc_file in modes.items():
    print(f"\n===== Processing mode: {mode} =====")
    print("Input:", fc_file)

    out_dir = Path(f"counts/from_raw_STAR_featureCounts_gencode_v37_{mode}")
    out_dir.mkdir(parents=True, exist_ok=True)

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

    print("Saved:", out_lnc)
    print("Shape:", count_out.shape)
    print("Library sizes:")
    print(count_out[sample_cols].sum(axis=0))

    paper_genes = [
        "AC024651.2", "AC026469.1", "AC129926.1", "AC244131.2",
        "AL031767.1", "AL732292.2", "LINC00324", "LINC02049",
        "ZNF436-AS1", "MIR17HG", "FP671120.7", "AC073140.2"
    ]

    print("\nSelected genes:")
    print(count_out[count_out["ID"].isin(paper_genes)].to_string(index=False))
