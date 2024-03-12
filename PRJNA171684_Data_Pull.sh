#!/bin/bash
#SBATCH --job-name=tissue_specific_Raw     # Job name
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=6                 # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=72:00:00                     # Time limit hrs:min:sec
#SBATCH --output=/scratch/tms51355/Taylor2024/rawdata/log.%j			# Location of standard output and error log files ##want an output in the work folder 
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALLd)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node



module load EDirect/21.4
module load SRA-Toolkit/3.0.3-gompi-2022a

bioproject="PRJNA171684"

OUTDIR="/scratch/tms51355/Taylor2024/rawdata/"                 

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

cd $OUTDIR

# Download SRA files using prefetch
prefetch "$bioproject"

sra_files=$(find "$OUTDIR"/ncbi/public/sra -name "*.sra")

# Convert SRA files to FASTQ format 
for sra_file in $sra_files; do
    fastq-dump --split-files "$sra_file"
done

echo "SRA files for BioProject $bioproject downloaded and converted to FASTQ."
echo tS