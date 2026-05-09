import pandas as pd
from pathlib import Path

res_file = Path("differential_expression_v37_MO_fraction/RUVSeq_edgeR_LRT_candidate_full_results.csv")
sig_file = Path("differential_expression_v37_MO_fraction/RUVSeq_edgeR_LRT_candidate_significant_paper_thresholds.csv")

res = pd.read_csv(res_file)
sig = pd.read_csv(sig_file)

paper_top10 = [
    "LINC02049", "AL031767.1", "AC244131.2", "LINC01954", "AL732292.2",
    "AC129926.1", "AC024651.2", "LINC02022", "AC026469.1", "LINC00324"
]

print("Number candidate significant:")
print(len(sig))

print("\nCandidate significant genes:")
if len(sig):
    print(sig[["ID", "logFC", "logCPM", "PValue", "FDR"]].to_string(index=False))
else:
    print("None")

subset = res[res["ID"].isin(paper_top10)].copy()

print("\nPaper top 10 genes in candidate full result:")
print(subset[["ID", "logFC", "logCPM", "PValue", "FDR"]].sort_values("ID").to_string(index=False))

print("\nPaper top 10 overlap with candidate significant:")
overlap = sorted(set(sig["ID"]) & set(paper_top10))
print(overlap)

print("\nPaper top 10 missing from candidate result:")
print(sorted(set(paper_top10) - set(res["ID"])))

# Also show top 30 by p-value
print("\nTop 30 candidate DE genes by PValue:")
print(res.sort_values("PValue")[["ID", "logFC", "logCPM", "PValue", "FDR"]].head(30).to_string(index=False))
