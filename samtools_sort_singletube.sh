#!/bin/bash
#SBATCH --job-name=TrimFilt  # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=28                  # Number of CPU cores per task
#SBATCH --mem=120gb                       # Job memory request
#SBATCH --time=72:00:00                     # Time limit hrs:min:sec
#SBATCH --output=trim_out         # Standard output log
#SBATCH --error=trim_error          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

 cd /scratch/tms51355/Taylor2024/singletube_NTS1

mkdir "Mapped_Data/hisat2_out"


for file in "Mapped_Data/demultiplexed/"*s.fastq*
do
	file2="${file:26:-9}"

if [ ! -f "Mapped_Data/hisat2_out/""$file2"".bam" ]; then

	module load SAMtools/1.16.1-GCC-11.3.0
	samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
	
fi
done

echo done