#!/bin/bash
#SBATCH --partition=nova
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=48
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH -J fc_bowtie1_v37
#SBATCH -o logs/featureCounts_bowtie1_v37_%j.out
#SBATCH -e logs/featureCounts_bowtie1_v37_%j.err
#SBATCH --hint=nomultithread

set -euo pipefail

cd /lustre/hdd/LAS/potoyan-lab/tuqa/project/paper_exact_pipeline

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /lustre/hdd/LAS/potoyan-lab/tuqa/conda_envs/rnaseq_env

mkdir -p counts/featureCounts_bowtie1_gencode_v37 logs

GTF="reference/gencode_v37/gencode.v37.annotation.gtf"
OUT="counts/featureCounts_bowtie1_gencode_v37/bowtie1_gencode_v37_featureCounts_s0_M_O_fraction.txt"

BAMS=(
  "bowtie1_alignments/C1_L-10_Aligned.sortedByCoord.out.bam"
  "bowtie1_alignments/C1_R+3_Aligned.sortedByCoord.out.bam"
  "bowtie1_alignments/C2_L-10_Aligned.sortedByCoord.out.bam"
  "bowtie1_alignments/C2_R+3_Aligned.sortedByCoord.out.bam"
  "bowtie1_alignments/C3_L-10_Aligned.sortedByCoord.out.bam"
  "bowtie1_alignments/C3_R+3_Aligned.sortedByCoord.out.bam"
)

echo "Running featureCounts on Bowtie1 BAMs"
echo "GTF: ${GTF}"
echo "Output: ${OUT}"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s 0 \
  -t exon \
  -g gene_id \
  -M \
  -O \
  --fraction \
  -a "${GTF}" \
  -o "${OUT}" \
  "${BAMS[@]}"

echo "featureCounts Bowtie1 completed."
