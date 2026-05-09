import pandas as pd
from pathlib import Path

m = pd.read_csv("metadata/sample_run_map.csv")

print("Number of rows in sample map:", len(m))

print("\nRuns per sample:")
print(m.groupby(["sample", "geo_sample"])["run"].count())

missing = []
for f in m["fastq"]:
    if not Path(f).exists():
        missing.append(f)

if missing:
    print("\nMissing FASTQ files:")
    for f in missing:
        print(f)
    raise SystemExit("ERROR: Some FASTQ files are missing.")
else:
    print("\nSUCCESS: All FASTQ files exist.")

print("\nFirst few rows:")
print(m.head())

