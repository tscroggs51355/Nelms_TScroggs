#!/bin/bash
#SBATCH --job-name=demulti
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=120gb
#SBATCH --time=72:00:00
#SBATCH --output=dm.%j.out
#SBATCH --error=dm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

source activate /home/tms51355/Fastq-Multx/
cd $SLURM_SUBMIT_DIR
module load fastp/0.23.2

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    # Skip if final output already exists
    if ls Mapped_Data/demultiplexed/${file2}_*_R2.fastq.gz 1>/dev/null 2>&1; then
        echo "Skipping ${file2} — already demultiplexed"
        continue
    fi

    (
        fastp -w 6 \
            -i "$file" \
            -I "Raw_Data/${file2}_R2_001.fastq.gz" \
            -o "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
            -O "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
            -A -Q -L -G \
            --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

        fastq-multx -b -B "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/CELSeq_barcodes" \
            -m 0 \
            "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
            "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
            -o "Mapped_Data/demultiplexed/${file2}_%_R2.fastq.gz" \
            "Mapped_Data/demultiplexed/${file2}_%.fastq.gz"
    ) &
done

wait  # allow all parallel jobs to finish
conda deactivate





























#!/bin/bash
#SBATCH --job-name=demulti
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=120gb
#SBATCH --time=72:00:00
#SBATCH --output=dm.%j.out
#SBATCH --error=dm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

source activate /home/tms51355/Fastq-Multx/
cd $SLURM_SUBMIT_DIR

module load fastp/0.23.2

logfile="processed_samples.log"
touch "$logfile"

for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

    if grep -q "$file2" "$logfile"; then
        echo "Already processed: $file2"
        continue
    fi

    if [ -f "Mapped_Data/demultiplexed/${file2}_R2.fastq.gz" ] && [ -f "Mapped_Data/demultiplexed/${file2}_R1.fastq.gz" ]; then
        echo "Skipping ${file2} — outputs exist"
        echo "$file2" >> "$logfile"
        continue
    fi

    (
        fastp -w 6 \
            -i "$file" \
            -I "Raw_Data/${file2}_R2_001.fastq.gz" \
            -o "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
            -O "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
            -A -Q -L -G \
            --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

        fastq-multx -b -B "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/CELSeq_barcodes" \
            -m 0 \
            "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
            "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
            -o "Mapped_Data/demultiplexed/${file2}_%_R2.fastq.gz" \
            "Mapped_Data/demultiplexed/${file2}_%.fastq.gz"

        echo "$file2" >> "$logfile"
    ) & 
done

wait  
conda deactivate