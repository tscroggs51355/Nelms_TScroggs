#!/bin/bash
#SBATCH --job-name=Mapping  # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=70gb                       # Job memory request
#SBATCH --time=48:00:00           # Time limit hrs:min:sec
#SBATCH --output=Mapping_output         # Standard output log
#SBATCH --error=Mapping_error          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node



module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0


for file in /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/Demultiplexed/*s.fastq*; do
    # Extract filename without path and remove .fastq or .fastq.gz
    filename=$(basename "$file")
    sample_name="${filename%.fastq.gz}"

    # Check if the *_unsorted.bam file already exists
    if [ -f "Mapped_Data/hisat2_out/${sample_name}_unsorted.bam" ]; then
        echo "Skipping ${sample_name} because *_unsorted.bam already exists."
        continue
    fi

    # Align with hisat2
    hisat2 -p 8 \
        -x /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/maize_tran \
        -U "$file" | samtools view -bS - > "Mapped_Data/hisat2_out/${sample_name}_unsorted.bam"

    # Sort the BAM file
    samtools sort -@ 8 "Mapped_Data/hisat2_out/${sample_name}_unsorted.bam" -o "Mapped_Data/hisat2_out/${sample_name}.bam"
done
