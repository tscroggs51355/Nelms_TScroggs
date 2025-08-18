#!/bin/bash
#SBATCH --job-name=UMICounts_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=UMICounts_NTS2_NTS2_1.%j.out
#SBATCH --error=UMICounts_NTS2_NTS2_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

cd /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2
mkdir -p Mapped_Data/bams
mkdir -p Mapped_Data/UMIcounts

for file in Mapped_Data/stringtie_out/*.bam; do
    file2=$(basename "$file" .bam)
    
    if [ ! -f "Mapped_Data/UMIcounts/${file2}.tsv" ]; then
        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 8 "$file" -o "Mapped_Data/bams/${file2}.sorted.bam"
        samtools index "Mapped_Data/bams/${file2}.sorted.bam"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
            -I "Mapped_Data/bams/${file2}.sorted.bam" \
            -S "Mapped_Data/UMIcounts/${file2}.tsv"
    fi
done




#!/bin/bash
#SBATCH --job-name=UMICounts_NTS2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=UMICounts_NTS2_%A_%a.out
#SBATCH --error=UMICounts_NTS2_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --array=0-$(($(wc -l < bam_list.txt)-1))
#SBATCH --export=NONE

# Load required modules
module load SAMtools/1.16.1-GCC-11.3.0
module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4

cd /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2
mkdir -p Mapped_Data/bams
mkdir -p Mapped_Data/UMIcounts

# Get the BAM file for this array index
bam_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" bam_list.txt)
file_base=$(basename "$bam_file" .bam)

# Only run if final output file is missing
if [ ! -f "Mapped_Data/UMIcounts/${file_base}.tsv" ]; then

    # Skip sorting if BAM index already exists
    if [ ! -f "Mapped_Data/bams/${file_base}.bai" ]; then
        samtools sort -@ 8 "$bam_file" -o "Mapped_Data/bams/${file_base}.bam"
        samtools index "Mapped_Data/bams/${file_base}.bam"
    fi

    umi_tools count \
        --per-gene \
        --gene-tag=XT \
        --assigned-status-tag=XS \
        -I "Mapped_Data/bams/${file_base}.bam" \
        -S "Mapped_Data/UMIcounts/${file_base}.tsv"

fi