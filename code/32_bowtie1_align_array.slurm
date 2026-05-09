#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --array=0-5
#SBATCH -J bowtie1_align
#SBATCH -o logs/bowtie1_align_%A_%a.out
#SBATCH -e logs/bowtie1_align_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p bowtie1_alignments logs

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
INDEX_PREFIX="reference/bowtie1_GRCh38/GRCh38_primary"

SAM_OUT="bowtie1_alignments/${SAMPLE}.sam"
BAM_OUT="bowtie1_alignments/${SAMPLE}_Aligned.sortedByCoord.out.bam"
LOG_OUT="bowtie1_alignments/${SAMPLE}_bowtie1.log"

echo "Sample: ${SAMPLE}"
echo "FASTQ: ${FASTQ}"
echo "Index: ${INDEX_PREFIX}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# Bowtie1 parameters:
# -S     output SAM
# -q     FASTQ input
# -p     threads
# -v 2   allow up to 2 mismatches in seed/alignment mode
# -m 20  suppress reads with >20 reportable alignments
# --best --strata report best-stratum alignments
#
# For gzipped FASTQ, feed with zcat.
zcat "${FASTQ}" | bowtie \
  -S \
  -q \
  -p "${SLURM_CPUS_PER_TASK}" \
  -v 2 \
  -m 20 \
  --best \
  --strata \
  "${INDEX_PREFIX}" \
  - \
  "${SAM_OUT}" \
  2> "${LOG_OUT}"

echo "Converting SAM to sorted BAM..."

samtools view -@ "${SLURM_CPUS_PER_TASK}" -bS "${SAM_OUT}" \
  | samtools sort -@ "${SLURM_CPUS_PER_TASK}" -o "${BAM_OUT}"

samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM_OUT}"

rm -f "${SAM_OUT}"

echo "Bowtie1 alignment completed for ${SAMPLE}"
ls -lh "bowtie1_alignments/${SAMPLE}"*
