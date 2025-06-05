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

for fastq_file in filtered/*.fastq.gz; do

    sample_name=$(basename "$fastq_file" .fastq.gz)

    output_file="filtered/${sample_name}_pUbUnspliced5.txt"

  zcat "$fastq_file" | sed -n '2~4p' | head -n 100000 | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTACGCCGCTCGTCCTCCCCCCC" | wc -l
done