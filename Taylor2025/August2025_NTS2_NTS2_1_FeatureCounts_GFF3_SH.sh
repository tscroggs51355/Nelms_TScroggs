#!/bin/bash
#SBATCH --job-name=UMICounts_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=12:00:00
#SBATCH --output=UMICounts_NTS2_NTS2_1.%j.out
#SBATCH --error=UMICounts_NTS2_NTS2_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE


module load GCC/11.3.0
module load ncurses/6.3-GCCcore-11.3.0
module load zlib/1.2.12-GCCcore-11.3.0
module load SAMtools/1.21-GCC-13.3.0
module load UMI-tools/1.1.4-foss-2023a

cd /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1


for file in "Mapped_Data/FeatureCount/"*.bam
do
	file2="${file:25:-18}"

	samtools sort -@ 8 "$file" -o "Mapped_Data/bams/""$file2"
	samtools index "Mapped_Data/bams/""$file2"
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "Mapped_Data/bams/""$file2" -S "Mapped_Data/UMIcounts/""${file2:0:-4}"".tsv"
done