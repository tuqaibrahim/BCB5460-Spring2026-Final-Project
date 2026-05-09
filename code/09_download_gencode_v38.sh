#!/usr/bin/env bash
set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

mkdir -p reference/gencode_v38
cd reference/gencode_v38

echo "Downloading GENCODE v38 GRCh38 primary assembly genome FASTA..."
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz

echo "Downloading GENCODE v38 primary assembly annotation GTF..."
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

echo "Unzipping reference files..."
gunzip -f GRCh38.primary_assembly.genome.fa.gz
gunzip -f gencode.v38.primary_assembly.annotation.gtf.gz

echo "Done."
ls -lh
