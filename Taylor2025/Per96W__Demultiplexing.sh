#!/bin/bash
#SBATCH --job-name=sampleD_run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80gb
#SBATCH --time=72:00:00
#SBATCH --output=sampleD.%j.out
#SBATCH --error=sampleD.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Activating conda environment..."
source activate /home/tms51355/Fastq-Multx/

cd $SLURM_SUBMIT_DIR
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Current directory: $SLURM_SUBMIT_DIR"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Loading fastp..."
module load fastp/0.23.2

file2="NTS2_NTS2_1_D_S39_L008"
R1="Raw_Data/${file2}_R1_001.fastq.gz"
R2="Raw_Data/${file2}_R2_001.fastq.gz"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input files:"
echo "    R1 = $R1"
echo "    R2 = $R2"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting fastp..."
start_fastp=$(date +%s)

fastp -w 6 \
    -i "$R1" \
    -I "$R2" \
    -o "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
    -O "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
    -A -Q -L -G \
    --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

end_fastp=$(date +%s)
runtime_fastp=$((end_fastp - start_fastp))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastp complete in ${runtime_fastp} seconds."

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output files from fastp:"
ls -lh Mapped_Data/demultiplexed/umi_${file2}_*.fastq.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting fastq-multx..."
start_multx=$(date +%s)

fastq-multx -b -B "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/CELSeq_barcodes" \
    -m 0 \
    "Mapped_Data/demultiplexed/umi_${file2}_R2.fastq.gz" \
    "Mapped_Data/demultiplexed/umi_${file2}_R1.fastq.gz" \
    -o "Mapped_Data/demultiplexed/${file2}_%_R2.fastq.gz" \
    "Mapped_Data/demultiplexed/${file2}_%.fastq.gz"

end_multx=$(date +%s)
runtime_multx=$((end_multx - start_multx))
echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastq-multx complete in ${runtime_multx} seconds."

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Demultiplexed files:"
ls -lh Mapped_Data/demultiplexed/${file2}_*.fastq.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Deactivating environment..."
conda deactivate

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample D pipeline finished."