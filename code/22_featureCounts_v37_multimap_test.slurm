#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --array=0-1
#SBATCH -J fc_v37_multi
#SBATCH -o logs/featureCounts_v37_multi_%A_%a.out
#SBATCH -e logs/featureCounts_v37_multi_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p counts/featureCounts_gencode_v37_multimap logs

GTF="reference/gencode_v37/gencode.v37.annotation.gtf"

BAMS=(
  "star_alignments/C1_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C1_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C2_R+3_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_L-10_Aligned.sortedByCoord.out.bam"
  "star_alignments/C3_R+3_Aligned.sortedByCoord.out.bam"
)

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    MODE="M"
    OUT="counts/featureCounts_gencode_v37_multimap/gencode_v37_featureCounts_s0_M.txt"
    EXTRA=(-M)
else
    MODE="M_fraction"
    OUT="counts/featureCounts_gencode_v37_multimap/gencode_v37_featureCounts_s0_M_fraction.txt"
    EXTRA=(-M --fraction)
fi

echo "Running featureCounts multimapper mode: ${MODE}"
echo "Output: ${OUT}"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s 0 \
  -t exon \
  -g gene_id \
  --primary \
  "${EXTRA[@]}" \
  -a "${GTF}" \
  -o "${OUT}" \
  "${BAMS[@]}"

echo "Done ${MODE}"
