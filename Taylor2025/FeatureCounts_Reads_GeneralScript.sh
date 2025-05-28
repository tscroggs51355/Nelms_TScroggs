#!/bin/bash
#SBATCH --job-name=featurecounts_March2025
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=FeatureCounts_March2025_stringtie2.%j.out
#SBATCH --error=FeatureCounts_March2025_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate subread-env

# Loop through each BAM file in the directory
for bam_file in /scratch/tms51355/Taylor2025/TS_March2025/Mapped_Data/hisat2_out/*.bam; do
    # Get the filename without extension
    filename=$(basename -- "$bam_file")
    filename_no_ext="${filename%.*}"

    # Run featureCounts for each BAM file
    featureCounts -T 4 -s 1 -a "/scratch/tms51355/Taylor2025/TS_March2025/Mapped_Data/stringtie_out/stringtie_merged.gtf" \
        -o "/scratch/tms51355/Taylor2025/TS_March2025/Mapped_Data/stringtie_out/${filename_no_ext}_read_counts.tab" \
        --readExtension5 500 -R BAM "$bam_file"
done