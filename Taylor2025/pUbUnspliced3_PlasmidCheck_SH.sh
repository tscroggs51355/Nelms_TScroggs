#!/bin/bash
#SBATCH --job-name=pUbUnspliced3
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=pUbUnspliced3.%j.out
#SBATCH --error=pUbUnspliced3.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

for fastq_file in filtered/*.fastq.gz; do

    sample_name=$(basename "$fastq_file" .fastq.gz)

    output_file="filtered/${sample_name}_pUbUnspliced3.txt"

    zcat "$fastq_file" | sed -n '2~4p' | grep "CCCTGTTGTTTGGTGTTACTTCTGCAGGTCGACTCTAGAGGATCCCCTCG" | \
        sed 's/CCCTGTTGTTTGGTGTTACTTCTGCAGGTCGACTCTAGAGGATCCCCTCG.*//' | \
        grep -E '^.{30,}$' | head -n 1000 | awk '{print substr($0, length($0) - 30 + 1)}' | \
        sort | uniq -c | sort -nr > "$output_file"
done

