#!/bin/bash
#SBATCH --job-name=featurecounts_March2025
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=FeatureCounts_June2025_stringtie2.%j.out
#SBATCH --error=FeatureCounts_June2025_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate subread-env

# Loop through each BAM file in the directory
for gtf_file in /scratch/tms51355/Taylor2025/TS_March2025/Mapped_Data/stringtie_out/*.gtf; do
    # Get the filename without extension
    filename=$(basename -- "$gtf_file")
    filename_no_ext="${filename%.*}"

# Run featureCounts on all BAM files together
featureCounts -T 6 -s 1 -a "Mapped_Data/stringtie_out/stringtie_merged.gtf" -o "Mapped_Data/stringtie_out/read_counts.tab" --readExtension5 500 -R BAM "Mapped_Data/hisat2_out/"*.bam