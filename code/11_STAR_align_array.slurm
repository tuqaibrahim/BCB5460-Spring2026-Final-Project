#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --array=0-5
#SBATCH -J STAR_align
#SBATCH -o logs/STAR_align_%A_%a.out
#SBATCH -e logs/STAR_align_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p star_alignments logs

SAMPLES=(
  "C1_L-10"
  "C1_R+3"
  "C2_L-10"
  "C2_R+3"
  "C3_L-10"
  "C3_R+3"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

FASTQ="trimmed_fastq/${SAMPLE}_cutadapt.fastq.gz"
GENOME_DIR="reference/STAR_gencode_v38"
OUT_PREFIX="star_alignments/${SAMPLE}_"

echo "Sample: ${SAMPLE}"
echo "FASTQ: ${FASTQ}"
echo "Genome index: ${GENOME_DIR}"
echo "Output prefix: ${OUT_PREFIX}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

STAR \
  --runThreadN "${SLURM_CPUS_PER_TASK}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${FASTQ}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --outSAMstrandField intronMotif \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04

echo "STAR alignment completed for ${SAMPLE}"

ls -lh star_alignments/${SAMPLE}_*
