#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=Merge_RC1November2024_stringtie2.%j.out
#SBATCH --error=Merge_RC1November2024_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

cd /scratch/tms51355/Taylor2024/RC1_Nov2024

ls -1 "Mapped_Data/stringtie_out/"*.gtf | gawk '{print $0}' > mergelist.txt

# Load StringTie module
module load StringTie/2.2.1-GCC-11.2.0

# Merge GTF files
stringtie --merge -p 6 -G /scratch/tms51355/Taylor2024/RC1_Nov2024/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o "Mapped_Data/stringtie_out/stringtie_merged.gtf" mergelist.txt