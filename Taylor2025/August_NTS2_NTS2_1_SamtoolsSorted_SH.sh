#!/bin/bash
#SBATCH --job-name=Mapping  # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=70gb                       # Job memory request
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --output=Mapping_output_Sam    # Standard output log
#SBATCH --error=Mapping_error_Sam     # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1

module load SAMtools/1.21-GCC-13.3.0

for file in "Mapped_Data/Demultiplexed/"*s.fastq.gz
do
	file2="${file:26:-9}"

	samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2""_sorted.bam"
	
done