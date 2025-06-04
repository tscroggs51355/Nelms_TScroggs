#!/bin/bash
#SBATCH --job-name=pUbSplice
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=pUbSplice.%j.out
#SBATCH --error=pUbSplice.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

for fastq_file in filtered/*.fastq.gz; do

    sample_name=$(basename "$fastq_file" .fastq.gz)

    output_file="filtered/${sample_name}_pUbSplice.txt"

    zcat "$fastq_file" | sed -n '2~4p' | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTCGACTCTAGAGGATCCCCTCG" | \
        sed 's/TCCACCCGTCGGCACCTCCGCTTCAAGGTCGACTCTAGAGGATCCCCTCG.*//' | \
        grep -E '^.{30,}$' | head -n 1000 | awk '{print substr($0, length($0) - 30 + 1)}' | \
        sort | uniq -c | sort -nr > "$output_file"
done