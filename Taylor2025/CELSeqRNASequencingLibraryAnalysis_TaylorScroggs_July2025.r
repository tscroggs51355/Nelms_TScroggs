# General Pipeline 
# R studio 
# CELSeq RNA Sequencing Library Analysis, Taylor Scroggs, 2025 

setwd("C:/Users/taylo/Desktop/June2025_Sequencing/")
annots = strsplit(read.table('Mapped_Data/stringtie_out/stringtie_merged.gtf', sep = '\t')[,9], '; ')
names(annots) = unlist(lapply(annots, function(xx) { xx[1] }))
names(annots) = sub('gene_id ', '', names(annots))
annots = annots[!duplicated(names(annots))]
annots = sub(';', '', sub(' ', '', unlist(lapply(annots, function(xx) { sub('.+ ', '', if (length(xx) == 3) { xx[3] } else { xx[1] }) }))))

files = dir('Mapped_Data/UMIcounts')
A0 = list()
for (f in files) {
  A0[[f]] = read.table(paste('Mapped_Data/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gn = unique(unlist(lapply(A0, rownames)))
A = matrix(NA, nrow = length(gn), ncol = length(files))
rownames(A) = gn
colnames(A) = files
for (f in files) {
  A[,f] = A0[[f]][match(gn,rownames(A0[[f]])),1]
}

## Need to change for each library as they might be named differently 
colnames(A) <- gsub("^Justin_|\\.bam\\.tsv$", "", colnames(A))
colnames(A) # just check what everything is 

A[is.na(A)] = 0
A = A[rowSums(A) > 0,]
rownames(A) = annots[rownames(A)]
dim(A)

A <- A[grepl("^Zm", rownames(A)), ]
dim(A)
 

for (g in unique(rownames(A)[duplicated(rownames(A))])) {
    i = which(rownames(A) == g)
    A[i[1],] = colSums(A[i,])
    A = A[-i[-1],]
}

dim(A)



####### UMI Counts and Genes 
# UMI Counts per Sample 
summary(colSums(A))

# Genes Per Sample 
summary(colSums(A > 0)) 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1132    9220   10903   10972   13886   15910


## Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/June2025_Sequencing/Mapped_Data/stringtie_out/read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- reads[, !grepl("_unsorted\\.bam$", colnames(reads))]
colnames(reads) <- sub("^Mapped_Data\\.hisat2_out\\.Justin_", "", colnames(reads))
colnames(reads) <- sub("\\.bam$", "", colnames(reads))

## Getting UMI Counts 
RperU = (reads[1,])/(colSums(A))
write.csv(RperU, "ReadsperUMICalculation_June2025Sequencing.csv")


## General Cluster Header 
#!/bin/bash
#SBATCH --job-name=pUbSplice
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=pUbSplice.%j.out
#SBATCH --error=pUbSplice.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

## Per sample, you will have four values 1) pUBSplice 2) pUbUnspliced3 3) pUbUnspliced5 4) InsertTag
## In this analysis, we want to see how often different parts of the plasmid sequence (spliced vs unspliced) were identified out of the TOTAL

## pUBSplice 
output_csv="demultiplexed/pUbSplice_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in demultiplexed/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 500000 | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTCGACTCTAGAGGATCCCCTCG" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done

## pUBUnspliced3 

output_csv="demultiplexed/pUbUnspliced3_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in demultiplexed/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 500000 | grep "CCCTGTTGTTTGGTGTTACTTCTGCAGGTCGACTCTAGAGGATCCCCTCG" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done

## pUBUnspliced5 

output_csv="demultiplexed/pUbUnspliced5_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in demultiplexed/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 500000 | grep "TCCACCCGTCGGCACCTCCGCTTCAAGGTACGCCGCTCGTCCTCCCCCCC" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done

## Insert Tag 

output_csv="demultiplexed/InsertTag_counts.csv"

echo "Sample,Count" > "$output_csv"

for fastq_file in demultiplexed/*.fastq.gz; do
    sample_name=$(basename "$fastq_file" .fastq.gz)

    count=$(zcat "$fastq_file" | sed -n '2~4p' | head -n 500000 | grep "GGGTGGGCGCG" | wc -l)

    echo "$sample_name,$count" >> "$output_csv"
done

## Data Wrangling for each of the four outputs above, need to change for each type of search 
filtered = "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data"

text_files_pUbSplice <- list.files(filtered, pattern = "_pUbSplice\\.txt$", full.names = TRUE)

data_list_pUbSplice <- lapply(text_files_pUbSplice, function(file) {
    sample_name <- gsub("_pUNos.txt", "", basename(file))

    # Check if file is empty
    if (file.info(file)$size > 0) {  
        df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence"))
        df$Sample <- sample_name
    } else {
        df <- data.frame(Count = NA, Sequence = NA, Sample = sample_name)  # Assign NA values
    }
    
    return(df)
})

final_data_pUbSplice <- do.call(rbind, data_list_pUbSplice)

# Save as CSV
write.csv(final_data_pUbSplice, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/pUbSplice_individual.csv", row.names = FALSE)

## In this analysis we are asking: 
## out of the subset of reads that have a sequence from the 3' plasmid transcript end
## what fraction have different specific insert sequences next to it? 

## Insert Tag Search 

for fastq_file in Mapped_Data/demultiplexed/*s.fastq.gz; do

    sample_name=$(basename "$fastq_file" .fastq.gz)

    zcat "$fastq_file" | sed -n '2~4p' | grep GGGTGGGCGCG | sed 's/GGGTGGGCGCG.*//' | \
        grep -E '^.{30,}$' | head -n 1000 | awk '{print substr($0, length($0) - 30 + 1)}' | \
        sort | uniq -c | sort -nr > "Mapped_Data/demultiplexed//${sample_name}.txt"


done

## Data Wrangling Insert Tag Search to make more readable 

library(tidyverse)

filtered = "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data"

text_files <- list.files("C:/Users/taylo/Desktop/TS_March2025/Mapped_Data", pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(text_files, function(file) {
    sample_name <- gsub(".txt", "", basename(file))  
    df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence")) 
    df$Sample <- sample_name 
    return(df)
})


final_data <- do.call(rbind, data_list)  

write.csv(final_data, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG_individual.csv", row.names = FALSE)


df_wide <- reshape(final_data, 
                   idvar = "Sequence", 
                   timevar = "Sample", 
                   direction = "wide")

colnames(df_wide) <- gsub("Count.", "", colnames(df_wide)) 

head(df_wide)

write.csv(df_wide, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG.csv", row.names = FALSE)


