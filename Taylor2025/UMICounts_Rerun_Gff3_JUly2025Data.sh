#!/bin/bash
#SBATCH --job-name=UMICounts_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=12:00:00
#SBATCH --output=UMICounts_NTS2_NTS2_1.%j.out
#SBATCH --error=UMICounts_NTS2_NTS2_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE


module load GCC/11.3.0
module load ncurses/6.3-GCCcore-11.3.0
module load zlib/1.2.12-GCCcore-11.3.0
module load SAMtools/1.21-GCC-13.3.0
module load UMI-tools/1.1.4-foss-2023a

cd /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/
 
for file in /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/FeatureCount/*.bam.featureCounts.bam; do
    # Get the base filename without the .bam.featureCounts.bam extension
    file2=$(basename "$file" .bam.featureCounts.bam)

    echo "Processing file: $file"

    # Check if the UMI count file already exists
    if [ ! -f "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/UMIcounts/${file2}.tsv" ]; then
        echo "UMI count file not found for $file2. Proceeding with sorting and counting."

        # Sort the BAM file using SAMtools
        echo "Sorting BAM file: $file"
        samtools sort -@ 8 "$file" -o "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/bams/${file2}.sorted.bam"
        samtools index "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/bams/${file2}.sorted.bam"
        echo "BAM file sorted and indexed: ${file2}.sorted.bam"

        # Count UMIs using umi_tools
        echo "Counting UMIs for $file2"
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
            -I "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/bams/${file2}.sorted.bam" \
            -S "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Mapped_Data/UMIcounts/${file2}.tsv"
        echo "UMI count file created: ${file2}.tsv"
    else
        echo "UMI count file already exists for $file2. Skipping..."
    fi
done

echo "Script completed."