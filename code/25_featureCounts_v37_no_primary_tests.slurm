#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --array=0-3
#SBATCH -J fc_v37_noPrimary
#SBATCH -o logs/featureCounts_v37_noPrimary_%A_%a.out
#SBATCH -e logs/featureCounts_v37_noPrimary_%A_%a.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p counts/featureCounts_gencode_v37_noPrimary logs

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
    MODE="M_noPrimary"
    OUT="counts/featureCounts_gencode_v37_noPrimary/gencode_v37_featureCounts_s0_M_noPrimary.txt"
    EXTRA=(-M)
elif [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    MODE="M_fraction_noPrimary"
    OUT="counts/featureCounts_gencode_v37_noPrimary/gencode_v37_featureCounts_s0_M_fraction_noPrimary.txt"
    EXTRA=(-M --fraction)
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
    MODE="M_O_noPrimary"
    OUT="counts/featureCounts_gencode_v37_noPrimary/gencode_v37_featureCounts_s0_M_O_noPrimary.txt"
    EXTRA=(-M -O)
else
    MODE="M_O_fraction_noPrimary"
    OUT="counts/featureCounts_gencode_v37_noPrimary/gencode_v37_featureCounts_s0_M_O_fraction_noPrimary.txt"
    EXTRA=(-M -O --fraction)
fi

echo "Running mode: ${MODE}"
echo "Output: ${OUT}"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s 0 \
  -t exon \
  -g gene_id \
  "${EXTRA[@]}" \
  -a "${GTF}" \
  -o "${OUT}" \
  "${BAMS[@]}"

echo "Done ${MODE}"

