#!/bin/bash
#SBATCH --job-name=combiningfiles_NTS2_NTS2_A_D
#SBATCH --partition=batch              # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=100gb                    # Job memory request
#SBATCH --time=72:00:00                     # Time limit hrs:min:sec
#SBATCH --output=combofiles_output
#SBATCH --error=combofiles_error
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=taylor.scroggs@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd /scratch/tms51355/Taylor2025/August_5_Sequencing_NTS2-NTS2-1/Raw_Data

cat NTS2-NTS2-1-A_S69_L008_R1_001.fastq.gz NTS2_NTS2_1_A_S36_L008_R1_001.fastq.gz > NTS2_NTS2_1_A_R1.fastq.gz

cat NTS2-NTS2-1-A_S69_L008_R2_001.fastq.gz NTS2_NTS2_1_A_S36_L008_R2_001.fastq.gz > NTS2_NTS2_1_A_R2.fastq.gz 

cat NTS2-NTS2-1-D_S70_L008_R1_001.fastq.gz NTS2_NTS2_1_D_S39_L008_R1_001.fastq.gz > NTS2_NTS2_1_D_R1.fastq.gz 

cat NTS2-NTS2-1-D_S70_L008_R2_001.fastq.gz NTS2_NTS2_1_D_S39_L008_R2_001.fastq.gz > NTS2_NTS2_1_D_R2.fastq.gz 