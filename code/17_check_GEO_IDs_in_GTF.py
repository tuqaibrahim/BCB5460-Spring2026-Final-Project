import pandas as pd
import re
from pathlib import Path

gtf_file = Path("reference/gencode_v38/gencode.v38.primary_assembly.annotation.gtf")
geo_file = Path("/lustre/hdd/LAS/potoyan-lab/tuqa/astronaut_exosomal_lncRNA_project/data/geo_processed/GSE193490_lncRNA_RuvSeq_EdgeR_counts.csv")

geo = pd.read_csv(geo_file)

ids_to_check = [
    # Paper top genes
    "LINC02049", "AL031767.1", "AC244131.2", "LINC01954", "AL732292.2",
    "AC129926.1", "AC024651.2", "LINC02022", "AC026469.1", "LINC00324",

    # Top GEO abundant genes
    "ZNF436-AS1", "MIR17HG", "MIRLET7BHG", "MIRLET7A1HG",
    "SNHG1", "MIR223HG", "LINC01783",

    # Accession-like IDs from GEO
    "FP671120.8", "FP236383.4", "AC073140.2"
]

# Also check the top 100 GEO IDs by total count
sample_cols = ["C1_L-10", "C1_R+3", "C2_L-10", "C2_R+3", "C3_L-10", "C3_R+3"]
geo["total"] = geo[sample_cols].sum(axis=1)
top_geo = geo.sort_values("total", ascending=False)["ID"].head(100).tolist()

ids_to_check = sorted(set(ids_to_check + top_geo))

def parse_attr(attr):
    out = {}
    for key, val in re.findall(r'(\S+) "([^"]+)"', attr):
        out[key] = val
    return out

records = []

print("Scanning GTF...")

with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, frame, attr = parts
        a = parse_attr(attr)

        row = {
            "feature": feature,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "gene_id": a.get("gene_id"),
            "gene_name": a.get("gene_name"),
            "gene_type": a.get("gene_type"),
            "transcript_id": a.get("transcript_id"),
            "transcript_name": a.get("transcript_name"),
            "transcript_type": a.get("transcript_type"),
        }

        values = {
            row["gene_id"], row["gene_name"],
            row["transcript_id"], row["transcript_name"]
        }

        hit_ids = [x for x in ids_to_check if x in values]

        for hit in hit_ids:
            r = row.copy()
            r["queried_ID"] = hit
            records.append(r)

hits = pd.DataFrame(records)

out = Path("counts/comparison_to_GEO/GEO_IDs_in_GENCODE_v38_GTF_hits.csv")
out.parent.mkdir(parents=True, exist_ok=True)
hits.to_csv(out, index=False)

print("\nNumber of queried IDs:", len(ids_to_check))
print("Number of hit rows:", len(hits))

if len(hits) > 0:
    print("\nHits grouped by queried ID and feature:")
    print(hits.groupby(["queried_ID", "feature"]).size().reset_index(name="n").to_string(index=False))

    print("\nFirst 100 hits:")
    print(hits.head(100).to_string(index=False))

missing = sorted(set(ids_to_check) - set(hits["queried_ID"])) if len(hits) else ids_to_check

print("\nQueried IDs not found in GTF gene_id/gene_name/transcript_id/transcript_name:")
print(missing)

print("\nSaved:")
print(out)
