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

output_csv="filtered/pUbSplice_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in filtered/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 100000 | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTCGACTCTAGAGGATCCCCTCG" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done

