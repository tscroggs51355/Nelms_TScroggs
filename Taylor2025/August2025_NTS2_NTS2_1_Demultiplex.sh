#!/bin/bash
#SBATCH --job-name=demulti
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80gb
#SBATCH --time=72:00:00
#SBATCH --output=dm.%j.out
#SBATCH --error=dm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE


module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/

echo "starting" 

# === Step 1: UMI handling and demultiplexing ===
for file in /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Raw_Data/*_R1*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2.fastq.gz//')
    
    echo "Processing file: $file2"

    if [ ! -f "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/umi_""$file2""_1s.fastq.gz" ]; then
        module load fastp/0.23.4-GCC-12.3.0

        echo "Running fastp on $file"
        fastp -w 6 -i "$file" -I "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Raw_Data/""$file2""_R2.fastq.gz" -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/umi_""$file2""_R1.fastq.gz" -O "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

        # Split read 2 file by CELseq barcodes. Require perfect match to barcode in expected location
        echo "Running fastq-multx on $file2"
        fastq-multx -b -B "CELSeq_barcodes" -m 0 "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/umi_""$file2""_R2.fastq.gz" "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/umi_""$file2""_R1.fastq.gz" "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Raw_Data/""$file2""_R2.fastq.gz" -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/""$file2""_%_R2.fastq.gz" "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/""$file2""_%.fastq.gz" "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/""$file2""_%_umi.fastq.gz"
    fi
done

conda deactivate

echo "moving to TRIM UMI containing reads for SRA upload"

# === Step 2: Trim UMI-containing reads for SRA upload ===
module load fastp/0.23.4-GCC-12.3.0

for file in /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/*s.fastq.gz; do
    file2=$(basename "$file" .fastq.gz)

    # Trim UMI containing read to only contain the UMI
    fastp -w 6 -B 10 -i "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/""$file2"".fastq.gz" -I "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/""$file2""_umi.fastq.gz" -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/SRA_upload/""$file2"".fastq.gz" -O "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/SRA_upload/""$file2""_umi.fastq.gz" -A -Q -L -G
done

echo "moving to polyA trimming" 

# === Step 3: PolyA trimming before alignment ===
for file in /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/*s.fastq.gz; do
    file2=$(basename "$file" .fastq.gz)

    if [ ! -f "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/hisat2_out/""$file2"".bam" ]; then
        fastp -w 6 -i "$file" -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA
    fi
done
