#!/bin/bash
#SBATCH --job-name=pUbUnspliced5
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=pUbUnspliced5.%j.out
#SBATCH --error=pUbUnspliced5.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

output_csv="demultiplexed/pUbUnspliced5_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in demultiplexed/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 500000 | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTACGCCGCTCGTCCTCCCCCCC" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done