#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=24G
#SBATCH --array=0-5
#SBATCH -J cutadapt_qc
#SBATCH -o logs/cutadapt_qc_%A_%a.out
#SBATCH -e logs/cutadapt_qc_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lusture/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

# Fix typo-safe package path if needed
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p trimmed_fastq cutadapt_logs qc/fastqc_trimmed logs

SAMPLES=(
  "C1_L-10"
  "C1_R+3"
  "C2_L-10"
  "C2_R+3"
  "C3_L-10"
  "C3_R+3"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

IN_FASTQ="raw_merged_fastq/${SAMPLE}.fastq.gz"
OUT_FASTQ="trimmed_fastq/${SAMPLE}_cutadapt.fastq.gz"
LOG_FILE="cutadapt_logs/${SAMPLE}_cutadapt.log"

echo "Sample: ${SAMPLE}"
echo "Input: ${IN_FASTQ}"
echo "Output: ${OUT_FASTQ}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

cutadapt \
  -j "${SLURM_CPUS_PER_TASK}" \
  -a AGATCGGAAGAGC \
  -a "A{10}" \
  -a "G{10}" \
  -m 15 \
  -o "${OUT_FASTQ}" \
  "${IN_FASTQ}" \
  > "${LOG_FILE}"

echo "Cutadapt finished for ${SAMPLE}"
echo "Running FastQC after trimming..."

fastqc \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir qc/fastqc_trimmed \
  "${OUT_FASTQ}"

echo "Done sample ${SAMPLE}"
