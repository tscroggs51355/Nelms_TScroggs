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
