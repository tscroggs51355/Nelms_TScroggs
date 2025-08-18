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

module load Miniforge3/24.11.3-0
MINIFORGE_BASE=$(conda info --base)
source ${MINIFORGE_BASE}/etc/profile.d/conda.sh
conda activate subread-env

mkdir "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/FeatureCount"

BAM_FILES=$(ls /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/hisat2_out/*s.bam | grep -Ev 'NTS2_NTS2_1_A|NTS2_NTS2_1_D')

featureCounts -T 8 -s 1 -a "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" \
  -t 'gene' -g 'ID' -o "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/FeatureCount/read_counts.tab" \
  --readExtension5 500 -R BAM $BAM_FILES