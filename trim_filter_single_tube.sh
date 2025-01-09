#!/bin/bash
#SBATCH --job-name=Samtools  # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=28                  # Number of CPU cores per task
#SBATCH --mem=120gb                       # Job memory request
#SBATCH --time=72:00:00                     # Time limit hrs:min:sec
#SBATCH --output=Samtools_sort_output         # Standard output log
#SBATCH --error=Samtools_sort_error          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd /scratch/tms51355/Taylor2024/singletube_NTS1

mkdir "Mapped_Data/hisat2_out"


for file in "Mapped_Data/demultiplexed/"*s.fastq*
do
	file2="${file:26:-9}"

if [ ! -f "Mapped_Data/hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.2.0
	fastp -w 4 -i "$file" -o "Mapped_Data/hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA
	
	module load HISAT2/2.2.1-gompi-2022aq4
	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 4 --dta -x /scratch/tms51355/Taylor2024/singletube_NTS1/maize_tran -U "Mapped_Data/hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"
	
	
	samtools sort -@ 4 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
	
	
fi
done
