#!/bin/bash
#SBATCH --job-name=Mapping_NTS2_NTS
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=MappingNTS2_1.%j.out
#SBATCH --error=MappingNTS2_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE



module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0


   for file in "Mapped_Data/Demultiplexed/NTS2_NTS2_1_D_34s.fastq.gz"
do
	file2="${file:26:-9}"

    # Align with hisat2
    hisat2 -p 8 \
        -x /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/maize_tran \
        -U "$file" | samtools view -bS - > "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"

    # Sort the BAM file
	samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"

done



###################################################################################
###################### FeautureCounts 


#!/bin/bash
#SBATCH --job-name=Feature_34s
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=FNTS2_1.%j.out
#SBATCH --error=FNTS2_1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

module load Miniforge3/24.11.3-0
MINIFORGE_BASE=$(conda info --base)
source ${MINIFORGE_BASE}/etc/profile.d/conda.sh
conda activate subread-env

BAM_FILE="/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/hisat2_out/NTS2_NTS2_1_D_34s.bam"

featureCounts -T 8 -s 1 -a "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" \
  -t 'gene' -g 'ID' -o "/scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Mapped_Data/FeatureCount/read_counts.tab" \
  --readExtension5 500 -R BAM "$BAM_FILE"






###################################################################################
### UMICounts 

#!/bin/bash
#SBATCH --job-name=UMICounts_NTS2_NTS2_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
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


for file in Mapped_Data/FeatureCount/NTS2_NTS2_1_D_*.bam
do
	file2="${file:25:-18}"

	samtools sort -@ 8 "$file" -o "Mapped_Data/bams/""$file2"
	samtools index "Mapped_Data/bams/""$file2"
	umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS -I "Mapped_Data/bams/""$file2" -S "Mapped_Data/UMIcounts/""${file2:0:-4}"".tsv"
done