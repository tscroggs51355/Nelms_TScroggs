#!/bin/bash
#SBATCH --job-name=stringtie2_out
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=singletube_stringtie2.%j.out
#SBATCH --error=singletube_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

cd /scratch/tms51355/Taylor2024/singletube_NTS1

mkdir "Mapped_Data/stringtie_out"

for file in "Mapped_Data/hisat2_out/"*.bam
do
        module load StringTie/2.2.1-GCC-11.2.0
        stringtie -p 4 -G /scratch/tms51355/Taylor2024/singletube_NTS1/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 --rf -o "Mapped_Data/stringtie_out/""${file:22:-4}"".gtf" "$file"
done