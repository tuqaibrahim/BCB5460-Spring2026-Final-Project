#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH -J STAR_index_v38
#SBATCH -o logs/STAR_index_v38_%j.out
#SBATCH -e logs/STAR_index_v38_%j.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p reference/STAR_gencode_v38 logs

FASTA="reference/gencode_v38/GRCh38.primary_assembly.genome.fa"
GTF="reference/gencode_v38/gencode.v38.primary_assembly.annotation.gtf"
INDEX_DIR="reference/STAR_gencode_v38"

echo "STAR version:"
STAR --version

echo "FASTA: $FASTA"
echo "GTF: $GTF"
echo "INDEX_DIR: $INDEX_DIR"
echo "CPUs: $SLURM_CPUS_PER_TASK"

STAR \
  --runThreadN "$SLURM_CPUS_PER_TASK" \
  --runMode genomeGenerate \
  --genomeDir "$INDEX_DIR" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 74

echo "STAR index build completed."
ls -lh "$INDEX_DIR"
