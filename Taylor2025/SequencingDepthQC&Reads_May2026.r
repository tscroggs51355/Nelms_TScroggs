## QC Looking at reads across NTS plates 
setwd("C:/Users/taylo/Desktop/TF Project Master Doucments")
NTS2 <- read.csv("NTS2_NTS2_1_Metadata.csv")
NTS3 <- read.csv("metadata_NTS3_NTS3_1.csv")

## Reads/UMI = NTS2-NTS2-1 
Reads_NTS2_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/ReadsperUMI_NTS2_NTS1_1.csv")
## Reads/UMI - NTS3-NTS3-1 
Reads_NTS3_NTS3_1 = read.csv("C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS3-NTS3-1/ReadsPerUMI_NTS3_NTS3_1.csv")

## Using the metadata from NTS2-NTS2 and NTS3_NTS3, need to change the colnames to match NTS2_NTS2_1_A_#s format 
well_to_sample_NTS2 <- setNames(NTS2$sample, NTS2$well)
cn <- colnames(Reads_NTS2_NTS2_1)

cn_new <- ifelse(cn %in% names(well_to_sample_NTS2),
                 well_to_sample_NTS2[cn],
                 cn)

colnames(Reads_NTS2_NTS2_1) <- cn_new


## 
well_to_sample_NTS3 <- setNames(NTS3$sample, NTS3$well)
cn1 <- colnames(Reads_NTS3_NTS3_1)

cn_2 <- ifelse(cn1 %in% names(well_to_sample_NTS3),
                 well_to_sample_NTS3[cn1],
                 cn1)

colnames(Reads_NTS3_NTS3_1) <- cn_2

## Then change NTS4-NTS4-1 and NTS1-NTS2-1 column names from  "NTS1.NTS2.1.D_96s" to Plate_Plate_Quadrant_CELSeqBarcode 

## Reads/UMI - NTS4-NTS4-1 
Reads_NTS4_NTS4_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/ReadsPerUMI_NTS4_NTS4_1.csv")
## Reads/UMI - NTS1-NTS2-1 
Reads_NTS1_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/ReadsPerUMI_NTS1_NTS2_1.csv")


Reads_NTS2_NTS2_1
Reads_NTS3_NTS3_1
Reads_NTS4_NTS4_1
Reads_NTS1_NTS2_1

extract_vals <- function(df) {
  as.numeric(df[ , -1])   # drop the first column ("X") and keep all 352 numeric wells
}

library(dplyr)
library(ggplot2)

plot_df <- bind_rows(
  data.frame(
    Plate = "NTS2_NTS2_1",
    RperU = extract_vals(Reads_NTS2_NTS2_1)
  ),
  data.frame(
    Plate = "NTS3_NTS3_1",
    RperU = extract_vals(Reads_NTS3_NTS3_1)
  ),
  data.frame(
    Plate = "NTS4_NTS4_1",
    RperU = extract_vals(Reads_NTS4_NTS4_1)
  ),

)

svg("C:/Users/taylo/Desktop/TF Project Master Doucments/ReadsperUMI_Boxplot_plates.svg", width = 8, height = 10)

ggplot(plot_df, aes(x = Plate, y = RperU, fill = Plate)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  theme_bw(base_size = 20) +   # increase overall font size
  labs(
    title = "Reads per UMI per Plate",
    y = "Reads per UMI"
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 24, face = "bold")
  )


dev.off()

# Percent Reads 
reads_2 = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/NTS4-NTS4-1-A_Seq/Mapped_Data/all_read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads_2) <- substr(colnames(reads_2), 72, nchar(colnames(reads_2)) - 4)
colnames(reads_2) <- sub("_S[0-9]+_", "_", colnames(reads_2))

suffixes <- sub(".*_", "", colnames(reads_2))
unique_tp <- unique(suffixes)
combined <- sapply(unique_tp, function(tp) {
  cols <- grep(paste0("_", tp, "$"), colnames(reads_2), value = TRUE)
  rowSums(reads_2[, cols, drop = FALSE])
})

colnames(combined) <- paste0("NTS4-NTS4-1-A_", unique_tp)
reads_2 <- as.data.frame(combined)


reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/November2025_Sequencing/Mapped_Data/all_read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 76, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads <- reads[, grep("^NTS4-NTS4-1", colnames(reads)), drop = FALSE]
reads_XX= cbind(reads,reads_2)

# 1. Subset only NTS4-NTS4-1 wells
nts4_cols <- grep("^NTS4-NTS4-1", colnames(reads_XX), value = TRUE)
reads_4 <- reads_XX[ , nts4_cols]

# 2. Total reads per well = Assigned + all Unassigned categories
total_reads_per_well <- colSums(reads_4)

# 3. Extract quadrant letter (A, B, C, D)
quadrant <- sub("^NTS4-NTS4-1-([A-D]).*", "\\1", names(total_reads_per_well))

# 4. Sum reads per quadrant
reads_per_quadrant <- tapply(total_reads_per_well, quadrant, sum)

# 5. Convert to percent of total
percent_quadrant <- 100 * reads_per_quadrant / sum(reads_per_quadrant)

percent_quadrant

## NTS2-NTS2-1 
setwd("C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August")


Augreads = read.table('C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August/Mapped_Data/Reads/all_read_counts.tab.summary',  header=T, sep = '\t', stringsAsFactors=F, row.names=1)
Julyreads = read.table('C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August/Mapped_Data/JulyReads/read_counts.tab.summary',  header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- cbind(Augreads,Julyreads)
write.csv(reads,"reads_updatedmapping.csv")

colnames(reads) <- substr(colnames(reads), 78, nchar(colnames(reads)) - 4)

colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))


colnames(reads)[1:192] <- sapply(colnames(reads)[1:192], function(x) substr(x, 9, nchar(x)))

reads_XX= reads

TFInfo = read.csv("TFInfo_Mod.csv")

metadata = read.csv("NTS2_NTS2_1_Metadata.csv")

reads_XX_clean <- reads_XX[ , !grepl("_(89|9[0-6])s$", names(reads_XX)) ]



# 1. Subset only NTS4-NTS4-1 wells
nts2_cols <- grep("^NTS2_NTS2_1", colnames(reads_XX_clean), value = TRUE)
reads_2 <- reads_XX_clean[ , nts2_cols]

# 2. Total reads per well = Assigned + all Unassigned categories
total_reads_per_well <- colSums(reads_2)

# 3. Extract quadrant letter (A, B, C, D)
quadrant <- sub("^NTS2_NTS2_1_([A-D]).*", "\\1", names(total_reads_per_well))

# 4. Sum reads per quadrant
reads_per_quadrant <- tapply(total_reads_per_well, quadrant, sum)

# 5. Convert to percent of total
percent_quadrant <- 100 * reads_per_quadrant / sum(reads_per_quadrant)

percent_quadrant

## NTS3_NTS3_1 
C:\Users\taylo\Desktop\2025_Nelms\August 2025\NTS3-NTS3-1\metadata_NTS3_NTS3_1.csv
metadata = read.csv("C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS3-NTS3-1/metadata_NTS3_NTS3_1.csv")
reads = read.table('C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS3-NTS3-1/reads/all_read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 70, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads_XX= reads
reads_XX_clean <- reads_XX[ , colnames(reads_XX) %in% metadata$sample ]

nts3_cols <- grep("^NTS3-NTS3-1", colnames(reads_XX_clean), value = TRUE)
reads_3 <- reads_XX_clean[ , nts3_cols]

# 2. Total reads per well = Assigned + all Unassigned categories
total_reads_per_well <- colSums(reads_3)

# 3. Extract quadrant letter (A, B, C, D)
quadrant <- sub("^NTS3-NTS3-1-([A-D])_.*", "\\1", names(total_reads_per_well))

# 4. Sum reads per quadrant
reads_per_quadrant <- tapply(total_reads_per_well, quadrant, sum)

# 5. Convert to percent of total
percent_quadrant <- 100 * reads_per_quadrant / sum(reads_per_quadrant)

percent_quadrant

## NTS1-NTS2-1 
reads = read.table('C:/Users/taylo/Desktop/TF Project Master Doucments/November2025_Sequencing/Mapped_Data/all_read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 76, nchar(colnames(reads)) - 4)
colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))
colnames(reads)<- gsub("\\.", "-", colnames(reads))
reads <- reads[, grep("^NTS1-NTS2-1", colnames(reads)), drop = FALSE]

nts_cols <- grep("^NTS1-NTS2-1", colnames(reads), value = TRUE)
reads_4 <- reads[ ,nts_cols]

# 2. Total reads per well = Assigned + all Unassigned categories
total_reads_per_well <- colSums(reads_4)

# 3. Extract quadrant letter (A, B, C, D)
quadrant <- sub("^NTS1-NTS2-1-([A-D])_.*", "\\1", names(total_reads_per_well))

# 4. Sum reads per quadrant
reads_per_quadrant <- tapply(total_reads_per_well, quadrant, sum)

# 5. Convert to percent of total
percent_quadrant <- 100 * reads_per_quadrant / sum(reads_per_quadrant)

percent_quadrant

### 5/26/2026 

# Updated regex: quadrant letter preceded by _ or . or -
# Handles: NTS2_NTS2_1_A_1s, NTS1.NTS2.1.A_9s, NTS3-NTS3-1-A_9s

datasets <- list(
  "NTS2_NTS2_1" = Reads_NTS2_NTS2_1,
  "NTS3_NTS3_1" = Reads_NTS3_NTS3_1,
  "NTS4_NTS4_1" = Reads_NTS4_NTS4_1,
  "NTS1_NTS2_1" = Reads_NTS1_NTS2_1
)

report <- do.call(rbind, lapply(names(datasets), function(ds_name) {
  df <- datasets[[ds_name]]
  df <- df[, colnames(df) != "X", drop = FALSE]
  
  # Match A/B/C/D preceded by _, ., or - and followed by _
  quadrants <- regmatches(colnames(df), regexpr("(?<=[_.\\-])[ABCD](?=_)", colnames(df), perl = TRUE))
  
  quad_means <- sapply(c("A", "B", "C", "D"), function(q) {
    cols <- colnames(df)[quadrants == q]
    if (length(cols) == 0) return(NA)
    mean(unlist(df[1, cols]), na.rm = TRUE)
  })
  
  data.frame(
    Dataset    = ds_name,
    Quadrant_A = round(quad_means["A"], 4),
    Quadrant_B = round(quad_means["B"], 4),
    Quadrant_C = round(quad_means["C"], 4),
    Quadrant_D = round(quad_means["D"], 4),
    row.names  = NULL
  )
}))

print(report)
write.csv(report, "quadrant_avg_reads_report.csv", row.names = FALSE)