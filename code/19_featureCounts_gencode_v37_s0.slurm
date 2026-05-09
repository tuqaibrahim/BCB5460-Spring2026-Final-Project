#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH -J featureCounts_v37
#SBATCH -o logs/featureCounts_v37_%j.out
#SBATCH -e logs/featureCounts_v37_%j.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

export CONDA_ENVS_PATH=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs
export CONDA_PKGS_DIRS=/lustre/hdd/LAS/potoyan-lab/tuqa/conda_pkgs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p counts/featureCounts_gencode_v37 logs

GTF="reference/gencode_v37/gencode.v37.annotation.gtf"
OUT="counts/featureCounts_gencode_v37/gencode_v37_featureCounts_s0.txt"

BAMS=(
  "star_alignments/C1_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C1_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_R+3_Aligned.sortedByCoord.out.bam"
)

echo "Running featureCounts with GENCODE v37, unstranded -s 0"
echo "GTF: ${GTF}"
echo "Output: ${OUT}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s 0 \
  -t exon \
  -g gene_id \
  --primary \
  -a "${GTF}" \
  -o "${OUT}" \
  "${BAMS[@]}"

echo "featureCounts v37 completed."
