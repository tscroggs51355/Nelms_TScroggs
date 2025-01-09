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

if [ ! -d "/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Data/hisat2_out" ]; then
  mkdir -p "/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Data/hisat2_out"
  echo "$(date): Created directory '/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Data/hisat2_out'."
fi

module load SAMtools/1.16.1-GCC-11.3.0

for file in "/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Data/demultiplexed/"*s.fastq*; do
  file2="${file:26:-9}"
  unsorted_bam="/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Data/hisat2_out/${file2}_unsorted.bam"
  sorted_bam="/scratch/tms51355/Taylor2024/singletube_NTS1/Mapped_Datata/hisat2_out/${file2}.bam"

  if [ ! -f "$sorted_bam" ]; then
    echo "$(date): Sorting file '$unsorted_bam'..."

    samtools sort -@ 8 "$unsorted_bam" -o "$sorted_bam"

    if [ $? -eq 0 ]; then
      echo "$(date): Sorting of '$unsorted_bam' completed successfully, output saved to '$sorted_bam'."
    else
      echo "$(date): Error sorting '$unsorted_bam'." >&2
    fi
  else
    echo "$(date): Sorted file '$sorted_bam' already exists, skipping..."
  fi
done

echo "$(date): Sorting process completed."
