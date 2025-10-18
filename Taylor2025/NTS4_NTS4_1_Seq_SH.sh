/work/bnlab/October_2025_Sequencing


nohup rsync /work/bnlab/October_2025_Sequencing /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq/Raw_Data


#!/bin/bash
#SBATCH --job-name=NTS4_NTS4_1_A_Demultiplexing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=NTS4_NTS4_1_A_Demultiplexing.%j.out
#SBATCH --error=NTS4_NTS4_1_A_Demultiplexing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
cd /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq/
module load fastp/0.23.4-GCC-13.2.0

for file in Raw_Data/*_R1_001.fastq.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//')

    echo "Processing $file2..."

    fastp -w 6 \
        -i "$file" \
        -I "Raw_Data/${file2}_R2_001.fastq.gz" \
        -o "Mapped_Data/Demultiplexed/umi_${file2}_R1.fastq.gz" \
        -O "Mapped_Data/Demultiplexed/umi_${file2}_R2.fastq.gz" \
        -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

    fastq-multx -b -B /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq/CELSeq_barcodes -m 1 \
        "Mapped_Data/Demultiplexed/umi_${file2}_R2.fastq.gz" \
        "Mapped_Data/Demultiplexed/umi_${file2}_R1.fastq.gz" \
        "Raw_Data/${file2}_R2_001.fastq.gz" \
        -o "Mapped_Data/Demultiplexed/${file2}_%_R2.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%_umi.fastq.gz"
done
#!/bin/bash
#SBATCH --job-name=NTS4_NTS4_1_A_Demultiplexing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=NTS4_NTS4_1_A_Demultiplexing.%j.out
#SBATCH --error=NTS4_NTS4_1_A_Demultiplexing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
cd /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq/Mapped_Data/
module load fastp/0.23.4-GCC-13.2.0

for file in Demultiplexed/*s.fastq.gz; do
        file2="${file:14:-9}"

    fastp -w 8 -B 10 -i "Demultiplexed/""$file2"".fastq.gz" -I "Demultiplexed/""$file2""_umi.fastq.gz" -o "SRA_upload/""$file2"".fastq.gz" -O "SRA_upload/""$file2""_umi.fastq.gz" -A -Q -L -G
done



for file in Demultiplexed/*s.fastq.gz; do
        file2="${file:14:-9}"

if [ ! -f "hisat2_out/""$file2"".bam" ]; then

        fastp -w 8 -i "$file" -o "hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

fi
done

### find /scratch/tms51355/Taylor2025/TS_March2025 -type f -exec touch {} +


#!/bin/bash
#SBATCH --job-name=NTS4_NTS4_1_A_Mapping_SH
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=NTS4_NTS4_1_A_Mapping.%j.out
#SBATCH --error=NTS4_NTS4_1_A_Mapping.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0
cd /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq

for file in "Mapped_Data/Demultiplexed/"*s.fastq*
do
        file2="${file:26:-9}"

         hisat2 -p 8 --dta -x /scratch/tms51355/Taylor2025/NTS4-NTS4-1-A-Seq/maize_tran -U "Mapped_Data/hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"
        samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
done
