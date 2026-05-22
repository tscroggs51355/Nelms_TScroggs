## NTS3-NTS3-1 Full Script on the cluster 

#####################################################################################################################################

#NTS3-NTS3-1 File Transfer 

#Location of files: 
/work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538

/work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538/NTS3* /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Raw_Data

nohup 
fpsync -n 8 -t $HOME/fpsync /work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538/NTS3* /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Raw_Data

fpsync -n 8 -t /home/tms51355/fpsync /work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538/NTS3* /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Raw_Data


fpsync -n 4 -t $HOME/fpsync /work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538/NTS3* /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Raw_Data
/home/tms51355/fpsync/parts/1756739629-2867270

rsync /work/bnlab/Aug_13_2025_Seq/24323-08-08192025_173538/NTS3* /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Raw_Data

#####################################################################################################################################

#!/bin/bash
#SBATCH --job-name=NTS3_NTS3_1_Demultiplexing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=12:00:00
#SBATCH --output=NTS3_NTS3_1_Demultiplexing.%j.out
#SBATCH --error=NTS3_NTS3_1_Demultiplexing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
cd /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/
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

    fastq-multx -b -B /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/CELSeq_barcodes -m 1 \
        "Mapped_Data/Demultiplexed/umi_${file2}_R2.fastq.gz" \
        "Mapped_Data/Demultiplexed/umi_${file2}_R1.fastq.gz" \
        "Raw_Data/${file2}_R2_001.fastq.gz" \
        -o "Mapped_Data/Demultiplexed/${file2}_%_R2.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%_umi.fastq.gz"
done

#!/bin/bash
#SBATCH --job-name=NTS3_NTS3_1_Demultiplexing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=12:00:00
#SBATCH --output=NTS3_NTS3_1_Demultiplexing.%j.out
#SBATCH --error=NTS3_NTS3_1_Demultiplexing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
cd /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Mapped_Data
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
###################################################################################
## Mapping ###

#!/bin/bash
#SBATCH --job-name=NTS3_NTS3_1_Mapping_SH
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=NTS3_NTS3_1_Mapping.%j.out
#SBATCH --error=NTS3_NTS3_1_Mapping.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0
cd /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq

for file in "Mapped_Data/Demultiplexed/"*s.fastq*
do
        file2="${file:26:-9}"

         hisat2 -p 8 --dta -x /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/maize_tran -U "Mapped_Data/hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"
        samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
done

###################################################################################
###################### FeautureCounts 


#!/bin/bash
#SBATCH --job-name=Feature_allreadcounts
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --output=NTS3_NTS3_1.%j.out
#SBATCH --error=NTS3_NTS3_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
MINIFORGE_BASE=$(conda info --base)
source ${MINIFORGE_BASE}/etc/profile.d/conda.sh
conda activate subread-env

cd /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq

mkdir "/scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Mapped_Data/FeatureCount"

BAM_FILES=$(ls /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Mapped_Data/hisat2_out/*s.bam)

featureCounts -T 8 -s 1 -a "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" \
  -t 'gene' -g 'ID' -o "/scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/Mapped_Data/FeatureCount/all_read_counts.tab" \
  --readExtension5 500 -R BAM $BAM_FILES


