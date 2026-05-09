#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --array=0-5
#SBATCH -J fastqc_raw
#SBATCH -o logs/fastqc_raw_%A_%a.out
#SBATCH -e logs/fastqc_raw_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p qc/fastqc_raw logs

SAMPLES=(
  "C1_L-10"
  "C1_R+3"
  "C2_L-10"
  "C2_R+3"
  "C3_L-10"
  "C3_R+3"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
FASTQ="raw_merged_fastq/${SAMPLE}.fastq.gz"

echo "Running FastQC for sample: ${SAMPLE}"
echo "FASTQ: ${FASTQ}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

fastqc \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --outdir qc/fastqc_raw \
  "${FASTQ}"

echo "Done FastQC for ${SAMPLE}"
