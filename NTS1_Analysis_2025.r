#NTS1 and single tube NTS1, files in different directories, but processing in parallel 

/scratch/tms51355/Taylor2024/NTS1 
/scratch/tms51355/Taylor2024/singletube_NTS1

## Moved files over on the xfer node in December 

### Have to move to these two directories: 
/scratch/tms51355/Taylor2024/NTS1 
/scratch/tms51355/Taylor2024/singletube_NTS1

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

## Changed line 48 for single tube samples ===  fastq-multx -b -B "/scratch/tms51355/Taylor2024/singletube_NTS1/CELSeq_barcodes" -m 0 "Mapped_Data/demultiplexed/umi_""$file2""_R2.fastq.gz" "Mapped_Data/demultiplexed/umi_""$file2""_R1.fastq.gz" -o "Mapped_Data/demultiplexed/""$file2""_%_R2.fastq.gz" "Mapped_Data/demultiplexed/""$file2""_%.fastq.gz"  for the single tube samples 


Step 2: Clean up files 

#### Just run in command line the following
###  -- have to be in the Mapped_Data/demultiplexed file path: 

Copy code
find "Mapped_Data/demultiplexed/" -name "umi_*" -print
find "Mapped_Data/demultiplexed/" -name "*s_R2*" -print


find "Mapped_Data/demultiplexed/" -name "umi_*" -delete
find "Mapped_Data/demultiplexed/" -name "*s_R2*" -delete


## Single Tube NTS1, Mapping Check 

fastp -w 4 -i Mapped_Data/demultiplexed/T7_S16_L002_1s.fastq.gz -o Mapped_Data/hisat2_out/T7_S16_L002_1s.fastq.gz -y -x -3 -a AAAAAAAAAAAA
fastp v0.23.2, time used: 7 seconds

The following have been reloaded with a version change:
  1) GCC/11.2.0 => GCC/11.3.0
  2) GCCcore/11.2.0 => GCCcore/11.3.0
  3) binutils/2.37-GCCcore-11.2.0 => binutils/2.38-GCCcore-11.3.0
  4) zlib/1.2.11-GCCcore-11.2.0 => zlib/1.2.12-GCCcore-11.3.0

1824749 reads; of these:
  1824749 (100.00%) were unpaired; of these:
    722662 (39.60%) aligned 0 times
    595879 (32.66%) aligned exactly 1 time
    506208 (27.74%) aligned >1 times
60.40% overall alignment rate
[bam_sort_core] merging from 0 files and 8 in-memory blocks...

The following have been reloaded with a version change:
  1) GCC/11.3.0 => GCC/11.2.0
  2) GCCcore/11.3.0 => GCCcore/11.2.0
  3) binutils/2.38-GCCcore-11.3.0 => binutils/2.37-GCCcore-11.2.0
  4) zlib/1.2.12-GCCcore-11.3.0 => zlib/1.2.11-GCCcore-11.2.0

Read1 before filtering:
total reads: 1467598
total bases: 221607298
Q20 bases: 203163330(91.6772%)
Q30 bases: 183739746(82.9123%)

Read1 after filtering:
total reads: 1438416
total bases: 166247464
Q20 bases: 161918389(97.396%)
Q30 bases: 154040242(92.6572%)

Filtering result:
reads passed filter: 1438416
reads failed due to low quality: 10570
reads failed due to too many N: 0
reads failed due to too short: 18380
reads failed due to low complexity: 232
reads with adapter trimmed: 948018
bases trimmed due to adapters: 49658567
reads with polyX in 3' end: 18089
bases trimmed in polyX tail: 520027

Duplication rate (may be overestimated since this is SE data): 24.9564%

JSON report: fastp.json
HTML report: fastp.html

fastp -w 4 -i Mapped_Data/demultiplexed/T8_S17_L002_1s.fastq.gz -o Mapped_Data/hisat2_out/T8_S17_L002_1s.fastq.gz -y -x -3 -a AAAAAAAAAAAA
fastp v0.23.2, time used: 6 seconds

The following have been reloaded with a version change:
  1) GCC/11.2.0 => GCC/11.3.0
  2) GCCcore/11.2.0 => GCCcore/11.3.0
  3) binutils/2.37-GCCcore-11.2.0 => binutils/2.38-GCCcore-11.3.0
  4) zlib/1.2.11-GCCcore-11.2.0 => zlib/1.2.12-GCCcore-11.3.0

1438416 reads; of these:
  1438416 (100.00%) were unpaired; of these:
    658411 (45.77%) aligned 0 times
    485117 (33.73%) aligned exactly 1 time
    294888 (20.50%) aligned >1 times
54.23% overall alignment rate



# Step 3: Stringtie 
#!/bin/bash
#SBATCH --job-name=stringtie2_out
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --output=stringtie2.%j.out
#SBATCH --error=stringtie2.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE


mkdir "Mapped_Data/stringtie_out"

for file in "Mapped_Data/hisat2_out/"*.bam
do
        module load StringTie/2.2.1-GCC-11.2.0
        stringtie -p 4 -G /scratch/tms51355/Taylor2024/August_2024_Sequencing_T/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 --rf -o "Mapped_Data/stringtie_out/""${file:22:-4}"".gtf" "$file"
done

# Step 4: Stringtie Merge 
#!/bin/bash
#SBATCH --job-name=stringtie2_merge
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --output=stringtie2merge.%j.out
#SBATCH --error=stringtie2merge.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR
ls -1 "Mapped_Data/stringtie_out/"*.gtf | gawk '{print $0}' > mergelist.txt

# Load StringTie module
module load StringTie/2.2.1-GCC-11.2.0

# Merge GTF files
stringtie --merge -p 4 -G /scratch/tms51355/Taylor2024/August_2024_Sequencing_T/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o "/scratch/tms51355/Taylor2024/August_2024_Sequencing_T/Mapped_Data/stringtie_out/stringtie_merge
d.gtf" mergelist.txt
