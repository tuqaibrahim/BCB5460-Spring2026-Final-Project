import pandas as pd
import subprocess
from pathlib import Path

m = pd.read_csv("metadata/sample_run_map.csv")

out_dir = Path("raw_merged_fastq")
out_dir.mkdir(exist_ok=True)

for sample, sub in m.groupby("sample", sort=False):
    out_fastq = out_dir / f"{sample}.fastq.gz"
    fastqs = sub["fastq"].tolist()

    print(f"\nMerging sample: {sample}")
    print(f"Output: {out_fastq}")
    for f in fastqs:
        print(f"  {f}")

    cmd = f"cat {' '.join(fastqs)} > {out_fastq}"
    subprocess.run(cmd, shell=True, check=True)

print("\nDone merging FASTQ files.")
