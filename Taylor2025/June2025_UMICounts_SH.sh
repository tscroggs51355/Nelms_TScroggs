#!/bin/bash
#SBATCH --job-name=UMICounts_June2025Sequencing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=UMICounts_June2025Seq_stringtie2.%j.out
#SBATCH --error=UMICounts_June2025Seq_stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE


cd /scratch/tms51355/Taylor2025/TS_June2025Sequencing
mkdir "Mapped_Data/bams"
mkdir "Mapped_Data/UMIcounts"

for file in "Mapped_Data/stringtie_out/"*.bam
do
    file2="${file:26:-18}"
    if [ ! -f "Mapped_Data/UMIcounts/${file2}.tsv" ]; then

        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort -@ 8 "$file" -o "Mapped_Data/bams/$file2"
        samtools index "Mapped_Data/bams/$file2"

        module load UMI-tools/1.1.2-foss-2022a-Python-3.10.4
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "Mapped_Data/bams/$file2" -S "Mapped_Data/UMIcounts/${file2}.tsv"
        # rm "$file"
    fi
done