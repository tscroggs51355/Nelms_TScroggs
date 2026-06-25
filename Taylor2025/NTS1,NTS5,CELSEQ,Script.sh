#Location of files: 
nohup wget -r -nH --cut-dirs=1 -P /work/bnlab/June2026Sequencing \
  --user=24323-Taylor.Scroggs --password=v8gRmM \
  ftp://38.122.175.98:2223/24323-15-05272026_151832

## Destination: 
/work/bnlab/June2026Sequencing/Batch_05292026_152726/deliver

nohup fpsync -n 8 -t $HOME/fpsync /work/bnlab/June2026Sequencing/Batch_05292026_152726/deliver /scratch/tms51355/Taylor2026/June2026Sequencing/Raw_Data


## Demultiplexing Step 1: 

#!/bin/bash
#SBATCH --job-name=PlateDemux
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=Plate_%A_%a.out
#SBATCH --error=Plate_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE
####
#!/bin/bash
#SBATCH --job-name=split
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --output=split.%j.out
#SBATCH --error=split.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

RAW_DIR=$(pwd)
RENAMED_DIR="${RAW_DIR}/renamed_fastqs"
mkdir -p "$RENAMED_DIR"

###############################################
# 1. Mapping table for 1T–24T → NTS names
###############################################

declare -A map

# Plate 1
map[1T]=NTS2_NTS3_1_A
map[2T]=NTS2_NTS3_1_B
map[3T]=NTS2_NTS3_1_C
map[4T]=NTS2_NTS3_1_D

# Plate 2
map[5T]=NTS1_NTS1_3_A
map[6T]=NTS1_NTS1_3_B
map[7T]=NTS1_NTS1_3_C
map[8T]=NTS1_NTS1_3_D

# Plate 3
map[9T]=NTS4_NTS5_1_A
map[10T]=NTS4_NTS5_1_B
map[11T]=NTS4_NTS5_1_C
map[12T]=NTS4_NTS5_1_D

# Plate 4
map[13T]=NTS3_NTS4_1_A
map[14T]=NTS3_NTS4_1_B
map[15T]=NTS3_NTS4_1_C
map[16T]=NTS3_NTS4_1_D

# Plate 5
map[17T]=NTS5_NTS5_1_A
map[18T]=NTS5_NTS5_1_B
map[19T]=NTS5_NTS5_1_C
map[20T]=NTS5_NTS5_1_D

# Plate 6
map[21T]=NTS5_NTS1_1_A
map[22T]=NTS5_NTS1_1_B
map[23T]=NTS5_NTS1_1_C
map[24T]=NTS5_NTS1_1_D

###############################################
# 2. Copy + rename FASTQs into renamed_fastqs/
###############################################

echo "Copying and renaming FASTQs…"

for f in *T_S*_L001_R*_001.fastq.gz; do
    prefix=$(echo "$f" | cut -d_ -f1)
    newprefix=${map[$prefix]}

    if [[ -z "$newprefix" ]]; then
        echo "Skipping $f — no mapping found"
        continue
    fi

    newname=$(echo "$f" | sed "s/^${prefix}/${newprefix}/")
    echo "$f → $newname"
    cp "$f" "${RENAMED_DIR}/${newname}"
done

###############################################
# 3. Split renamed files into plate directories
###############################################

cd "$RENAMED_DIR"

mkdir -p Plate_{1..6}

mv NTS2_NTS3_1_* Plate_NTS2_NTS3_1/
mv NTS1_NTS1_3_* Plate_NTS1_NTS1_3/
mv NTS4_NTS5_1_* Plate_NTS4_NTS5_1/
mv NTS3_NTS4_1_* Plate_NTS3_NTS4_1/
mv NTS5_NTS5_1_* Plate_NTS5_NTS5_1/
mv NTS5_NTS1_1_* Plate_NTS5_NTS1_1/

echo "Files split into plate directories."

#!/bin/bash
#SBATCH --job-name=PlateDemux
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=12:00:00
#SBATCH --output=Plate_%A_%a.out
#SBATCH --error=Plate_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

###############################################
# 1. Load modules and activate environment
###############################################
module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
module load fastp/0.23.4-GCC-13.2.0

###############################################
# 2. Determine which plate this array task runs
###############################################
BASE_DIR=/scratch/tms51355/Taylor2026/June2026Sequencing/Raw_Data/renamed_fastqs
PLATE_DIR=${BASE_DIR}/Plate_${SLURM_ARRAY_TASK_ID}

echo "Processing plate directory: $PLATE_DIR"

# Safety check
if [ ! -d "$PLATE_DIR" ]; then
    echo "ERROR: Plate directory $PLATE_DIR does not exist"
    exit 1
fi

cd "$PLATE_DIR"

###############################################
# 3. Create output directories inside each plate
###############################################
mkdir -p Mapped_Data/Demultiplexed
mkdir -p Mapped_Data/SRA_upload
mkdir -p Mapped_Data/hisat2_out
mkdir -p Mapped_Data/FeatureCount

###############################################
# 4. Run your existing workflow unchanged
###############################################
for file in *_R1_001.fastq.gz; do
    filename=$(basename "$file")
    file2=$(echo "$filename" | sed 's/_R1.*//')

    echo "Processing $file2..."

    fastp -w 6 \
        -i "$file" \
        -I "${file2}_R2_001.fastq.gz" \
        -o "Mapped_Data/Demultiplexed/umi_${file2}_R1.fastq.gz" \
        -O "Mapped_Data/Demultiplexed/umi_${file2}_R2.fastq.gz" \
        -A -Q -L -G --umi --umi_loc read2 --umi_len 10 --umi_prefix UMI

    fastq-multx -b -B /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/CELSeq_barcodes -m 1 \
        "Mapped_Data/Demultiplexed/umi_${file2}_R2.fastq.gz" \
        "Mapped_Data/Demultiplexed/umi_${file2}_R1.fastq.gz" \
        "${file2}_R2_001.fastq.gz" \
        -o "Mapped_Data/Demultiplexed/${file2}_%_R2.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%.fastq.gz" \
           "Mapped_Data/Demultiplexed/${file2}_%_umi.fastq.gz"
done

echo "Plate ${SLURM_ARRAY_TASK_ID} complete."

#### 
#!/bin/bash
#SBATCH --job-name=PostDemux
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --output=PostDemux_%A_%a.out
#SBATCH --error=PostDemux_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=taylor.scroggs@uga.edu
#SBATCH --export=NONE

###############################################
# 1. Load modules and activate environment
###############################################
module load Miniforge3/24.11.3-0
source activate /home/tms51355/Fastq-Multx/
module load fastp/0.23.4-GCC-13.2.0

###############################################
# 2. Determine which plate this array task runs
###############################################
BASE_DIR=/scratch/tms51355/Taylor2026/June2026Sequencing/Raw_Data/renamed_fastqs
PLATE_DIR=${BASE_DIR}/Plate_${SLURM_ARRAY_TASK_ID}/Mapped_Data

echo "Processing: $PLATE_DIR"

if [ ! -d "$PLATE_DIR" ]; then
    echo "ERROR: $PLATE_DIR does not exist"
    exit 1
fi

cd "$PLATE_DIR"

###############################################
# 3. Make sure output directories exist
###############################################
mkdir -p SRA_upload
mkdir -p hisat2_out

###############################################
# 4. Loop 1 — SRA_upload fastp processing
###############################################
echo "Running SRA_upload fastp..."

for file in Demultiplexed/*s.fastq.gz; do
    # Extract sample name from path
    file2="${file:14:-9}"

    fastp -w 8 -B 10 \
        -i "Demultiplexed/${file2}.fastq.gz" \
        -I "Demultiplexed/${file2}_umi.fastq.gz" \
        -o "SRA_upload/${file2}.fastq.gz" \
        -O "SRA_upload/${file2}_umi.fastq.gz" \
        -A -Q -L -G
done

###############################################
# 5. Loop 2 — hisat2_out polyA trimming
###############################################
echo "Running hisat2_out trimming..."

for file in Demultiplexed/*s.fastq.gz; do
    file2="${file:14:-9}"

    if [ ! -f "hisat2_out/${file2}.bam" ]; then
        fastp -w 8 \
            -i "$file" \
            -o "hisat2_out/${file2}.fastq.gz" \
            -y -x -3 -a AAAAAAAAAAAA
    fi
done

echo "Plate ${SLURM_ARRAY_TASK_ID} complete."

## 
#!/bin/bash
#SBATCH --job-name=AlignCount
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=12:00:00
#SBATCH --output=AlignCount_%A_%a.out
#SBATCH --error=AlignCount_%A_%a.err
#SBATCH --export=NONE

###############################################
# 1. Load modules
###############################################
module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.21-GCC-13.3.0
module load Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate subread-env

###############################################
# 2. Plate directory
###############################################
BASE=/scratch/tms51355/Taylor2026/June2026Sequencing/Raw_Data/renamed_fastqs
PLATE_DIR=${BASE}/Plate_${SLURM_ARRAY_TASK_ID}/Mapped_Data

cd "$PLATE_DIR"

###############################################
# 3. HISAT2 alignment
###############################################
for file in Demultiplexed/*s.fastq.gz; do
    fname=$(basename "$file")
    sample=${fname:0:-9}

    hisat2 -p 8 --dta \
        -x /scratch/tms51355/Taylor2025/NTS3_NTS3_1_Seq/maize_tran \
        -U "hisat2_out/${sample}.fastq.gz" \
        | samtools view -bS - > "hisat2_out/${sample}_unsorted.bam"

    samtools sort -@ 8 \
        "hisat2_out/${sample}_unsorted.bam" \
        -o "hisat2_out/${sample}.bam"

    rm "hisat2_out/${sample}_unsorted.bam"
done

###############################################
# 4. featureCounts
###############################################
BAMS=$(ls hisat2_out/*s.bam)

featureCounts -T 8 -s 1 \
  -a "/scratch/tms51355/Taylor2025/July2025Sequencing_NTS2/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" \
  -t gene -g ID \
  -o "FeatureCount/Plate_${SLURM_ARRAY_TASK_ID}_counts.tab" \
  --readExtension5 500 -R BAM $BAMS

echo "Plate ${SLURM_ARRAY_TASK_ID} complete."

## Plate_1 Analysis 


setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_1")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_1/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS2_NTS3_1_A_", bc), paste0("NTS2_NTS3_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS2_NTS3_1_C_", bc), paste0("NTS2_NTS3_1_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS2_NTS3_1_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS2-NTS3", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS2-NTS3"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS2_NTS3_1.csv")

## Normalized Data in Z, filtred data in Z_filted for NTS2-NTS3-1 
A2 <- read.csv("A2_NTS2_NTS3_1_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_1/Plate_1_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS2_NTS3_1.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS2_NTS3_1.csv", row.names = FALSE)

summary(colSums(A2[,grep("A", colnames(A2))]))
summary(colSums(A2[,grep("B", colnames(A2))]))
summary(colSums(A2[,grep("C", colnames(A2))]))
summary(colSums(A2[,grep("D", colnames(A2))]))
summary(colSums(A2))

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    350   10209   14654   17164   21132   78403 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    117    8593   12260   14616   18402   48875 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1928    7646   12472   14584   19121   60100 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     41    7386   12481   16048   20302   65264 
> summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     41    8070   12906   15603   19955   78403 


## Plate_2 Analysis, NTS1-NTS1-3 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_2")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_2/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS1_NTS1_3_A_", bc), paste0("NTS1_NTS1_3_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS1_NTS1_3_C_", bc), paste0("NTS1_NTS1_3_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS1_NTS1_3_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS1-NTS1", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS1-NTS1"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS1_NTS1_3.csv")

## Normalized Data in Z, filtred data in Z_filted for NTS2-NTS3-1 
A2 <- read.csv("A2_NTS1_NTS1_3_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_2/Plate_2_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS1_NTS1_3.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS1_NTS1_3.csv", row.names = FALSE)

summary(colSums(A2[,grep("A", colnames(A2))]))
summary(colSums(A2[,grep("B", colnames(A2))]))
summary(colSums(A2[,grep("C", colnames(A2))]))
summary(colSums(A2[,grep("D", colnames(A2))]))
summary(colSums(A2))

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   3276    7166   10642   13157   15751   54029 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1589    6489   10320   13162   16408   50913 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   2710    9414   15782   22050   27107   95305 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1835    7402   11740   16424   20740   72295 
> summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1589    7409   11500   16198   19369   95305 

## Plate_3 Analysis, NTS4-NTS5-1 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_3")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_3/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS4_NTS5_1_A_", bc), paste0("NTS4_NTS5_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS4_NTS5_1_C_", bc), paste0("NTS4_NTS5_1_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS4_NTS5_1_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS4-NTS5", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS4-NTS5"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS4_NTS5_1.csv")

A2 <- read.csv("A2_NTS4_NTS5_1_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_3/Plate_3_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS4_NTS5_1.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS4_NTS5_1.csv", row.names = FALSE)

summary(colSums(A2[,grep("A", colnames(A2))]))
summary(colSums(A2[,grep("B", colnames(A2))]))
summary(colSums(A2[,grep("C", colnames(A2))]))
summary(colSums(A2[,grep("D", colnames(A2))]))
summary(colSums(A2))

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    849    4281    7994    8537   11978   21290 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   2768    8314   14438   16472   21058   53936 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1539    6384   11692   12894   18698   51305 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1040    5548   10598   13123   17831   43366 
> summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    849    5789   10879   12756   17741   53936

## Plate_4 Analysis, NTS3-NTS4-1 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_4")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_4/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS3_NTS4_1_A_", bc), paste0("NTS3_NTS4_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS3_NTS4_1_C_", bc), paste0("NTS3_NTS4_1_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS3_NTS4_1_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS3-NTS4", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS3-NTS4"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS3_NTS4_1.csv")

A2 <- read.csv("A2_NTS3_NTS4_1_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_4/Plate_4_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS3_NTS4_1.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS3_NTS4_1.csv", row.names = FALSE)

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    397    4994    9451   13405   19176   67262 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    806    7755   12528   15739   18711   48863 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    224    8013   13212   15753   20539   59303 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1522    8471   16962   18822   23393   53303 

### Plate_5, NTS5-NTS5-1

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_5")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_5/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS5_NTS5_1_A_", bc), paste0("NTS5_NTS5_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS5_NTS5_1_C_", bc), paste0("NTS5_NTS5_1_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS5_NTS5_1_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS5-NTS5", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS5-NTS5"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS5_NTS5_1.csv")

A2 <- read.csv("A2_NTS5_NTS5_1_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_5/Plate_5_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS5_NTS5_1.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS5_NTS5_1.csv", row.names = FALSE)

summary(colSums(A2[,grep("A", colnames(A2))]))
summary(colSums(A2[,grep("B", colnames(A2))]))
summary(colSums(A2[,grep("C", colnames(A2))]))
summary(colSums(A2[,grep("D", colnames(A2))]))
summary(colSums(A2))

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   3314   10529   16618   22792   31961  170376 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   4038   10631   17243   21894   29863  104514 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   5115   15516   25898   30167   39329   87730 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   4312   13542   23800   26346   37623   68814 
> summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   3314   11970   19918   25300   34622  170376 

## Plate_6, NTS5-NTS1-1 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_6")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv")

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

files = dir('UMIcounts')
files <- files[!grepl("_[0-8]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_6/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files
colnames(A2) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2))
colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}

A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(9:96, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(9:96, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS5_NTS1_1_A_", bc), paste0("NTS5_NTS1_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS5_NTS1_1_C_", bc), paste0("NTS5_NTS1_1_D_", bc))
}))


Wellorder <- c(AB, CD)
Wellorder <- Wellorder[Wellorder %in% colnames(A2)]
A2 <- A2[, Wellorder, drop = FALSE]



metaID <- data.frame(
  sample = colnames(A2),
  well   = meta$well)

  sample_lookup <- setNames(metadata$SampleID, metadata$Array.Well)
metaID$SampleID <- sample_lookup[metaID$well]


sample_to_well <- setNames(metaID$well, metaID$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
write.csv(A2, "A2_NTS5_NTS1_1_PreNormalization.csv")

metadata = read.csv("Array/metadata_TFAtlas.csv")
metadata <- metadata[grepl("NTS5-NTS1", metadata$Assay.Plate), ]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS5-NTS1"]]

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(4)) > 0, ]
write.csv(Z, "Z_NTS5_NTS1_1.csv")

A2 <- read.csv("A2_NTS5_NTS1_1_PreNormalization.csv",
               row.names = 1,
               check.names = FALSE)
well_to_sample <- setNames(metaID$sample, metaID$well)
colnames(A2) <- well_to_sample[colnames(A2)]

reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_6/Plate_6_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 12, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads


valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)

svg(" ReadsPerUMIBarplot.svg", width = 8, height = 8)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7)
dev.off()

write.csv(RperU, "ReadsPerUMI_NTS5_NTS1_1.csv")

qc_df <- data.frame(
  well            = colnames(A2),
  UMIcounts       = colSums(A2),
  gene_counts     = colSums(A2 > 0),
  ReadsperUMI     = as.numeric(RperU[colnames(A2)]),
  UMI_per_gene    = colSums(A2) / colSums(A2 > 0),
  row.names       = NULL
)

write.csv(qc_df, "QC_summary_NTS5_NTS1_1.csv", row.names = FALSE)

summary(colSums(A2[,grep("A", colnames(A2))]))
summary(colSums(A2[,grep("B", colnames(A2))]))
summary(colSums(A2[,grep("C", colnames(A2))]))
summary(colSums(A2[,grep("D", colnames(A2))]))
summary(colSums(A2))

> summary(colSums(A2[,grep("A", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    350   10209   14654   17164   21132   78403 
> summary(colSums(A2[,grep("B", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    117    8593   12260   14616   18402   48875 
> summary(colSums(A2[,grep("C", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1928    7646   12472   14584   19121   60100 
> summary(colSums(A2[,grep("D", colnames(A2))]))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     41    7386   12481   16048   20302   65264 
> summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     41    8070   12906   15603   19955   78403

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk/Plate_6")

QC_NTS5_NTS1_1 = read.csv("QC_summary_NTS5_NTS1_1.csv")

library(dplyr)

QC_NTS5_NTS1_1 %>%
  mutate(quadrant = sub(".*_([A-D])_.*", "\\1", well)) %>%
  group_by(quadrant) %>%
  summarise(
    mean_ReadsperUMI = mean(ReadsperUMI),
    mean_UMIcounts   = mean(UMIcounts),
    mean_gene_counts = mean(gene_counts),
    mean_UMI_per_gene = mean(UMI_per_gene)
  )

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk")

QC_NTS5_NTS1_1 = read.csv("QC_summary_NTS5_NTS1_1.csv")
QC_NTS2_NTS3_1 = read.csv("QC_summary_NTS2_NTS3_1.csv")
QC_NTS5_NTS5_1 = read.csv("QC_summary_NTS5_NTS5_1.csv")
QC_NTS3_NTS4_1 = read.csv("QC_summary_NTS3_NTS4_1.csv")
QC_NTS3_NTS4_1$reads_per_gene <- NULL
QC_NTS1_NTS1_3 = read.csv("QC_summary_NTS1_NTS1_3.csv")
QC_NTS4_NTS5_1 = read.csv("QC_summary_NTS4_NTS5_1.csv")

QC_ALL <- rbind(QC_NTS1_NTS1_3,QC_NTS2_NTS3_1,QC_NTS3_NTS4_1,QC_NTS4_NTS5_1,QC_NTS5_NTS1_1,QC_NTS5_NTS5_1)

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/NTS1,NTS5Bulk") 

metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Metadata_FullPlate_Jan2026.csv")

### 1. Extract quadrant from SampleName
metadata$Quadrant <- sub(".*-([A-D])_.*", "\\1", metadata$SampleName)

### 2. Build metadata lookup
meta_lookup <- metadata[, c("Assay.Plate", "Array.Well", "SampleID", "SampleName", "Quadrant")]

A2_NTS2_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS2_NTS2_1.csv", row.names = 1, check.names = FALSE)
A2_NTS3_NTS3_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS3_NTS3_1.csv", row.names = 1, check.names = FALSE)
A2_NTS4_NTS4_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS4_NTS4_1_PreNormalization.csv", row.names = 1, check.names = FALSE)
A2_NTS1_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS1_NTS2_1_PreNormalization.csv", row.names = 1, check.names = FALSE)

A2_list <- list(
  NTS2_NTS2_1 = A2_NTS2_NTS2_1,
  NTS3_NTS3_1 = A2_NTS3_NTS3_1,
  NTS4_NTS4_1 = A2_NTS4_NTS4_1,
  NTS1_NTS2_1 = A2_NTS1_NTS2_1
)

files <- list.files(
  pattern = "^ReadsPerUMI_.*\\.csv$",
  full.names = TRUE
)

reads_list <- lapply(files, read.csv)

# Name each list element by the file name (without .csv)
names(reads_list) <- gsub("\\.csv$", "", basename(files))
## Helper: fix column names for one ReadsPerUMI table using metadata
fix_RperU_cols <- function(df, plate_label_meta) {

  # 1) Normalize column names: dots → dashes
  sample_names <- gsub("\\.", "-", colnames(df))

  # 2) Get metadata for this plate (Assay.Plate)
  plate_meta <- metadata[metadata$Assay.Plate == plate_label_meta, ]

  # 3) Build mapping: SampleName -> Array.Well
  map <- setNames(plate_meta$Array.Well, plate_meta$SampleName)

  # 4) Translate SampleName (from file) to wells
  new_cols <- unname(map[sample_names])

  # 5) Apply new column names
  colnames(df) <- new_cols
  df
}

## NTS4-NTS4-1
reads_list$ReadsPerUMI_NTS4_NTS4_1 <-
  fix_RperU_cols(reads_list$ReadsPerUMI_NTS4_NTS4_1, "NTS4-NTS4")

## NTS1-NTS2-1
reads_list$ReadsPerUMI_NTS1_NTS2_1 <-
  fix_RperU_cols(reads_list$ReadsPerUMI_NTS1_NTS2_1, "NTS1-NTS2")

reads_list$ReadsPerUMI_NTS2_NTS2_1 <- reads_list$ReadsPerUMI_V5_NTS2_NTS2_1
reads_list$ReadsPerUMI_V5_NTS2_NTS2_1 <- NULL

make_QC <- function(A2, RperU_vec) {

  data.frame(
    well         = colnames(A2),
    UMIcounts    = colSums(A2),
    gene_counts  = colSums(A2 > 0),
    ReadsperUMI  = as.numeric(RperU_vec[colnames(A2)]),
    UMI_per_gene = colSums(A2) / colSums(A2 > 0),
    stringsAsFactors = FALSE
  )
}

QC_list <- list()

for (plate in names(A2_list)) {

  Rname <- paste0("ReadsPerUMI_", plate)

  # Extract numeric ReadsPerUMI row
  Rvec <- as.numeric(reads_list[[Rname]][1, ])
  names(Rvec) <- colnames(reads_list[[Rname]])

  # Build QC table
  QC_list[[plate]] <- make_QC(
    A2 = A2_list[[plate]],
    RperU_vec = Rvec
  )
}

library(dplyr)

QC_FIVE <- bind_rows(QC_list, .id = "Plate")

QC_FIVE$SampleName <- paste(QC_FIVE$Plate, QC_FIVE$well, sep = "_")
QC_FIVE <- QC_FIVE[, c("SampleName", setdiff(names(QC_FIVE), c("SampleName", "Plate", "well")))]

QC_Combined = rbind(QC_FIVE, QC_ALL)

# 1. Rename 'well' to 'SampleName' in QC_ALL
names(QC_ALL)[names(QC_ALL) == "well"] <- "SampleName"

# 2. Reorder QC_ALL columns to match QC_FIVE
QC_ALL <- QC_ALL[, names(QC_FIVE)]

# 3. Combine
QC_COMBINED <- rbind(QC_FIVE, QC_ALL)

write.csv(QC_COMBINED, "QC_June2026_TFProject.csv")