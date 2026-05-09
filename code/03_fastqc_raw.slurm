#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH -J fastqc_raw
#SBATCH -o logs/fastqc_raw_%j.out
#SBATCH -e logs/fastqc_raw_%j.err

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p qc/fastqc_raw

fastqc \
  --threads 6 \
  --outdir qc/fastqc_raw \
  raw_merged_fastq/*.fastq.gz

echo "FastQC raw completed."
ls -lh qc/fastqc_raw
