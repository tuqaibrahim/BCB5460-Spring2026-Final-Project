#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH -J bowtie1_index
#SBATCH -o logs/bowtie1_index_%j.out
#SBATCH -e logs/bowtie1_index_%j.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p reference/bowtie1_GRCh38 logs

FASTA="reference/gencode_v38/GRCh38.primary_assembly.genome.fa"
INDEX_PREFIX="reference/bowtie1_GRCh38/GRCh38_primary"

echo "Bowtie version:"
bowtie --version

echo "Building Bowtie1 index..."
echo "FASTA: ${FASTA}"
echo "Index prefix: ${INDEX_PREFIX}"

bowtie-build --threads "${SLURM_CPUS_PER_TASK}" "${FASTA}" "${INDEX_PREFIX}"

echo "Bowtie1 index completed."
ls -lh reference/bowtie1_GRCh38
