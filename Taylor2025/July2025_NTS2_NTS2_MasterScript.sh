## Demultiplexing 
# Activate your Conda environment here if needed
source activate /home/tms51355/Fastq-Multx/

cd $SLURM_SUBMIT_DIR


for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

if [ ! -f "Mapped_Data/demultiplexed/""$file2""_dT-1s.fastq.gz" ]; then
module load fastp/0.23.2
	fastp -w 4 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -O "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

	fastq-multx -b -B "/scratch/tms51355/Taylor2025/TS_March2025/CELSeq_barcodes" -m 0 "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -o "Mapped_Data/demultiplexed/""$file2""_%_R2.fastq.gz" "Mapped_Data/demultiplexed/""$file2""_%.fastq.gz"  

fi
done
conda deactivate

#############################################################################333
### Mapping 
cd /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2

mkdir "Mapped_Data/hisat2_out"


for file in "Mapped_Data/demultiplexed/"*s.fastq*
do
        file2="${file:26:-9}"

if [ ! -f "Mapped_Data/hisat2_out/""$file2"".bam" ]; then

        module load fastp/0.23.2-GCC-11.2.0
        fastp -w 8 -i "$file" -o "Mapped_Data/hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA

        module load HISAT2/2.2.1-gompi-2022a
        module load SAMtools/1.16.1-GCC-11.3.0
        hisat2 -p 8  -x /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/maize_tran -U "Mapped_Data/hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
fi
done






###################################################################################
###################### FeautureCounts 
#!/bin/bash
#SBATCH --job-name=F_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=12:00:00
#SBATCH --output=FNTS2_NTS2_1.%j.out
#SBATCH --error=FNTS2_NTS2_1.%j.err
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






###################################################################################
### UMICounts 





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