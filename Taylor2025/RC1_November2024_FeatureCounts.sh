#!/bin/bash
#SBATCH --job-name=featurecounts_RC1_November2024
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=FeatureCounts_RC1_November2024_stringtie2.%j.out
#SBATCH --error=FeatureCounts_RC1_November2024_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

cd /scratch/tms51355/Taylor2024/RC1_Nov2024

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate subread-env

featureCounts -T 6 -s 1 \
  -a "Mapped_Data/stringtie_out/stringtie_merged.gtf" \
  -o "Mapped_Data/stringtie_out/read_counts.tab" \
  --readExtension5 500 -R BAM Mapped_Data/hisat2_out/*.bam