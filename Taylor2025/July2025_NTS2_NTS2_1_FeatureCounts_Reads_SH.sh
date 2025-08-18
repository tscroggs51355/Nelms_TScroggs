#!/bin/bash
#SBATCH --job-name=featurecounts_July2025_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=FeatureCounts_July2025_NTS2_NTS2_1_stringtie2.%j.out
#SBATCH --error=FeatureCounts_July2025_NTS2_NTS2_1_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate subread-env

# Run featureCounts once for all BAMs
featureCounts -T 8 -s 1 \
  -a /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/stringtie_out/stringtie_merged.gtf \
  -o /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/stringtie_out/all_samples_read_counts.tab \
  --readExtension5 500 -R BAM \
  /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/hisat2_out/*.bam