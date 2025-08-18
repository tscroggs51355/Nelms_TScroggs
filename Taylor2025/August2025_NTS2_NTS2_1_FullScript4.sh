## August 2025 Full Script 

###################################################################################
## Demultiplexing Samples 

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

###################################################################################
## Mapping ###

module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0


   for file in "Mapped_Data/Demultiplexed/"*s.fastq.gz
do
	file2="${file:26:-9}"

    # Check if the *_unsorted.bam file already exists
    if [ -f "Mapped_Data/hisat2_out/${sample_name}_unsorted.bam" ]; then
        echo "Skipping ${sample_name} because *_unsorted.bam already exists."
        continue
    fi

    # Align with hisat2
    hisat2 -p 8 \
        -x /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/maize_tran \
        -U "$file" | samtools view -bS - > "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"

    # Sort the BAM file
	samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2""_sorted.bam"

done

###################################################################################
###################### FeautureCounts 
module load Miniforge3/24.11.3-0
MINIFORGE_BASE=$(conda info --base)
source ${MINIFORGE_BASE}/etc/profile.d/conda.sh
conda activate subread-env

mkdir "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/FeatureCount"

BAM_FILES=$(ls /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/hisat2_out/*s_sorted.bam)

featureCounts -T 8 -s 1 -a "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" \
  -t 'gene' -g 'ID' -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/FeatureCount/read_counts.tab" \
  --readExtension5 500 -R BAM $BAM_FILES






###################################################################################
### UMICounts 

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

cd /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1


for file in "Mapped_Data/FeatureCount/"*.bam
do
	file2="${file:25:-18}"

	samtools sort -@ 8 "$file" -o "Mapped_Data/bams/""$file2"
	samtools index "Mapped_Data/bams/""$file2"
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "Mapped_Data/bams/""$file2" -S "Mapped_Data/UMIcounts/""${file2:0:-4}"".tsv"
done