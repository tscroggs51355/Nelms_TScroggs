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

module load SAMtools/1.16.1-GCC-11.3.0
for file in Mapped_Data/hisat2_out/*_unsorted.bam
do
    base_name="${file%_unsorted.bam}"
    
    samtools sort -@ 8 "$file" -o "${base_name}.bam"
    
done