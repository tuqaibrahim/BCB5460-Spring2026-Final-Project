#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --array=0-2
#SBATCH -J featureCounts_strand
#SBATCH -o logs/featureCounts_strand_%A_%a.out
#SBATCH -e logs/featureCounts_strand_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p counts/featureCounts_strand_test logs

STRANDS=(0 1 2)
S=${STRANDS[$SLURM_ARRAY_TASK_ID]}

GTF="reference/gencode_v38/gencode.v38.primary_assembly.annotation.gtf"
OUT="counts/featureCounts_strand_test/gencode_v38_featureCounts_s${S}.txt"

BAMS=(
  "star_alignments/C1_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C1_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_R+3_Aligned.sortedByCoord.out.bam"
)

echo "Running featureCounts with strandedness -s ${S}"
echo "GTF: ${GTF}"
echo "Output: ${OUT}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s "${S}" \
  -t exon \
  -g gene_id \
  --primary \
  -a "${GTF}" \
  -o "${OUT}" \
  "${BAMS[@]}"

echo "featureCounts completed for -s ${S}"
