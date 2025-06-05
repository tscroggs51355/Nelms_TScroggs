#!/bin/bash
#SBATCH --job-name=GGGTGGGCGCG_NTS1_NTS1_1_PlasmidCheck
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=PlasmidCheck_NTS1_NTS1_1.%j.out
#SBATCH --error=PlasmidCheck_NTS1_NTS1_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

for fastq_file in Mapped_Data/demultiplexed/*s.fastq.gz; do

    sample_name=$(basename "$fastq_file" .fastq.gz)

    zcat "$fastq_file" | sed -n '2~4p' | grep GGGTGGGCGCG | sed 's/GGGTGGGCGCG.*//' | \
        grep -E '^.{30,}$' | head -n 1000 | awk '{print substr($0, length($0) - 30 + 1)}' | \
        sort | uniq -c | sort -nr > "Mapped_Data/demultiplexed//${sample_name}.txt"


done
