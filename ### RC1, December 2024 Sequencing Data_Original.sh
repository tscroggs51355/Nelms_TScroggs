### RC1, December 2024 

Copy Files over while on the Xfer Node -- from work/bnlab to scratch 

cp CELSeq_barcodes /scratch/tms51355/Taylor2024/RC1_Nov2024/
cp maize*  /scratch/tms51355/Taylor2024/RC1_Nov2024/
cp Zm*  /scratch/tms51355/Taylor2024/RC1_Nov2024/

mkdir Raw_Data
mkdir Mapped_Data
mkdir Mapped_Data/demultiplexed

Step 1: ## Demultiplex Samples 

#!/bin/bash
#SBATCH --job-name=demulti      # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=4                  # Number of CPU cores per task
#SBATCH --mem=80gb                          # Job memory request
#SBATCH --time=5:00:00                     # Time limit hrs:min:sec
#SBATCH --output=dm.%j.out         # Standard output log
#SBATCH --error=dm.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALLd)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

# Activate your Conda environment here if needed
source activate /home/tms51355/Fastq-Multx/
cd $SLURM_SUBMIT_DIR


for file in Raw_Data/*_R1_*.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//' | sed 's/_R2_001.fastq.gz//')

if [ ! -f "Mapped_Data/demultiplexed/""$file2""_dT-1s.fastq.gz" ]; then
module load fastp/0.23.2
	fastp -w 4 -i "$file" -I "Raw_Data/""$file2""_R2_001.fastq.gz" -o "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -O "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

	fastq-multx -b -B "/scratch/tms51355/Taylor2024/NTS1/CELSeq_barcodes" -m 0 "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -o "Mapped_Data/demultiplexed/""$file2""_%_R2.fastq.gz" "Mapped_Data/demultiplexed/""$file2""_%.fastq.gz"  

fi
done
conda deactivate

Step 2: Clean up files 

#### Just run in command line the following
###  -- have to be in the Mapped_Data/demultiplexed file path: 

Copy code
find "Mapped_Data/demultiplexed/" -name "umi_*" -print
find "Mapped_Data/demultiplexed/" -name "*s_R2*" -print


find "Mapped_Data/demultiplexed/" -name "umi_*" -delete
find "Mapped_Data/demultiplexed/" -name "*s_R2*" -delete

Step  3: Trim and Filter 

#!/bin/bash
#SBATCH --job-name=TrimFilt  # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=28                  # Number of CPU cores per task
#SBATCH --mem=120gb                       # Job memory request
#SBATCH --time=72:00:00                     # Time limit hrs:min:sec
#SBATCH --output=trim_out         # Standard output log
#SBATCH --error=trim_error          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd $SLURM_SUBMIT_DIR

mkdir "Mapped_Data/hisat2_out"


for file in "Mapped_Data/demultiplexed/"*s.fastq*
do
	file2="${file:26:-9}"

if [ ! -f "Mapped_Data/hisat2_out/""$file2"".bam" ]; then

	module load fastp/0.23.2-GCC-11.2.0
	fastp -w 4 -i "$file" -o "Mapped_Data/hisat2_out/""$file2"".fastq.gz" -y -x -3 -a AAAAAAAAAAAA
	
	module load HISAT2/2.2.1-gompi-2022a
	module load SAMtools/1.16.1-GCC-11.3.0
	hisat2 -p 4 --dta -x /scratch/tms51355/Taylor2024/RC1_Nov2024/maize_tran -U "Mapped_Data/hisat2_out/""$file2"".fastq.gz" | samtools view -bS -> "Mapped_Data/hisat2_out/""$file2""_unsorted.bam"
	
	module load SAMtools/1.16.1-GCC-11.3.0
	samtools sort -@ 8 "Mapped_Data/hisat2_out/""$file2""_unsorted.bam" -o "Mapped_Data/hisat2_out/""$file2"".bam"
	
fi
done



Step 4: Clean up files 
	find "Mapped_Data/hisat2_out/" -name "*fastq.gz" -delete
	find "Mapped_Data/hisat2_out/" -name "*_unsorted.bam" -delete


File Check: 

High multimapping, on par with October 2023

fastp -w 4 -i Mapped_Data/demultiplexed/RC1_S22_L002_10s.fastq.gz -o Mapped_Data/hisat2_out/RC1_S22_L002_10s.fastq.gz -y -x -3 -a AAAAAAAAAAAA
fastp v0.23.2, time used: 5 seconds

The following have been reloaded with a version change:
  1) GCC/11.2.0 => GCC/11.3.0
  2) GCCcore/11.2.0 => GCCcore/11.3.0
  3) binutils/2.37-GCCcore-11.2.0 => binutils/2.38-GCCcore-11.3.0
  4) zlib/1.2.11-GCCcore-11.2.0 => zlib/1.2.12-GCCcore-11.3.0

1657820 reads; of these:
  1657820 (100.00%) were unpaired; of these:
    256313 (15.46%) aligned 0 times
    154901 (9.34%) aligned exactly 1 time
    1246606 (75.20%) aligned >1 times
84.54% overall alignment rate
[bam_sort_core] merging from 0 files and 8 in-memory blocks...

The following have been reloaded with a version change:
  1) GCC/11.3.0 => GCC/11.2.0
  2) GCCcore/11.3.0 => GCCcore/11.2.0
  3) binutils/2.38-GCCcore-11.3.0 => binutils/2.37-GCCcore-11.2.0
  4) zlib/1.2.12-GCCcore-11.3.0 => zlib/1.2.11-GCCcore-11.2.0

Read1 before filtering:
total reads: 1036302
total bases: 156481602
Q20 bases: 142183897(90.863%)
Q30 bases: 127685078(81.5975%)

Read1 after filtering:
total reads: 994381
total bases: 118530501
Q20 bases: 115175085(97.1692%)
Q30 bases: 109645033(92.5036%)

Filtering result:
reads passed filter: 994381
reads failed due to low quality: 9957
reads failed due to too many N: 0
reads failed due to too short: 31415
reads failed due to low complexity: 549
reads with adapter trimmed: 528457
bases trimmed due to adapters: 33312231
reads with polyX in 3' end: 19317
bases trimmed in polyX tail: 782805

Duplication rate (may be overestimated since this is SE data): 38.4783%

JSON report: fastp.json
HTML report: fastp.html

fastp -w 4 -i Mapped_Data/demultiplexed/RC1_S22_L002_11s.fastq.gz -o Mapped_Data/hisat2_out/RC1_S22_L002_11s.fastq.gz -y -x -3 -a AAAAAAAAAAAA
fastp v0.23.2, time used: 4 seconds

The following have been reloaded with a version change:
  1) GCC/11.2.0 => GCC/11.3.0
  2) GCCcore/11.2.0 => GCCcore/11.3.0
  3) binutils/2.37-GCCcore-11.2.0 => binutils/2.38-GCCcore-11.3.0
  4) zlib/1.2.11-GCCcore-11.2.0 => zlib/1.2.12-GCCcore-11.3.0

994381 reads; of these:
  994381 (100.00%) were unpaired; of these:
    182346 (18.34%) aligned 0 times
    122930 (12.36%) aligned exactly 1 time
    689105 (69.30%) aligned >1 times
81.66% overall alignment rate
[bam_sort_core] merging from 0 files and 8 in-memory blocks...


gunzip RC1_S22_L002_1s.fastq.gz | sed -n '2~4p' | head -n 100000 | sort | uniq -c | sort -nr | head -n 10

   1650 GGCGGTGGATACCTAGGCACCCAGAGACGAGGAAGGGCGTAGCAAGCGACGAAATGCTTCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAAC
   1308 GGGGGCTCGAAGACGATCAGATACCGTCCTAGTCTCAACCATAAACGATGCCGACCAGGGATCAGCGGGTGTTACTAATAGGACCCCGCTGGCACCTTATGAGAAATCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAA 
   ## 18s Ribosomal 
   1116 GGGGGCATTCGTATTTCATAGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTCTCAACCATA
   1046 TGCGGTGGATACCTAGGCACCCAGAGACGAGGAAGGGCGTAGCAAGCGACGAAATGCTTCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAAC
    610 GGGCGTAGCAAGCGACGAAATGCTTCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAACCTGGCGAACTGAAACATCTTAGTAGCCAGAGGAA
    565 GCCCTATCAACTTTCGATGGTAGGATAGGGGCCTACCATGGTGGTGACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGG
	## 18s Ribosomal 
    470 GCCCTGTCAACTTTCGATGGTAGGATAGGGGCCTACCATGGTGGTGACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGG
	# 18s Ribosomal 
    436 TTGCGGTGGATACCTAGGCACCCAGAGACGAGGAAGGGCGTAGCAAGCGACGAAATGCTTCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAA
    434 TCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAACCTGGCGAACTGAAACATCTTAGTAGCCAGAGGAAAAAAAAAAAAAAAAAAAAAAAAAA
    423 GGGCGGTGGATACCTAGGCACCCAGAGACGAGGAAGGGCGTAGCAAGCGACGAAATGCTTCGGGGAGTTGAAAATAAGCATAGATCCGGAGATTCCCAAATAGGTCAACCTTTTGAACTGCCTGCTGAATCCATGAGCAGGCAAGAGACAA