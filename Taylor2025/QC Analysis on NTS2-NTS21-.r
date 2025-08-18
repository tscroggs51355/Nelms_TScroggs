## QC Analysis: 
load("UMI count data.rda")
"UMI count data.rda"
####
setwd("C:/Users/taylo/Desktop/2025_Nelms/July2025/NTS2_Seq")


Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/2025_Nelms/July2025/NTS2_Seq/Mapped_Data/stringtie_out/all_samples_read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)
colnames(reads) <- substr(colnames(reads), 78, nchar(colnames(reads)))
print(colnames(reads))
colnames(reads) <- gsub("_S\\d+_L008|\\.bam$", "", colnames(reads))

colnames(reads) = metadata$well[match(colnames(reads), metadata$sample)]

## Getting UMI Counts 
RperU = (reads[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
colnames(RperU) <- colnames(reads)
RperU_matrix <- as.matrix(RperU)
head(RperU_matrix)
RperU_matrix <- RperU_matrix[, match(metadata$well, colnames(RperU_matrix))]
head(RperU_matrix)

#olnames(RperU_matrix) <- metadata$well[match(colnames(RperU_matrix), metadata$sample)]

RperU_df <- as.data.frame(RperU_matrix)
head(RperU_df)
RperU_df$Sample <- rownames(RperU_df)
head(RperU_df)

ConvertMatrixtoPlate <- function(xx) {
  # We need exactly 352 values to fill the plate (16 rows x 22 columns)
  # Create a matrix with 16 rows and 22 columns first (without the first and last columns)
  mat_data <- matrix(xx, nrow = 16, ncol = 22, byrow = TRUE)
  
  # Now create the full 16x24 plate with 1st and last columns as NA
  mat <- matrix(NA, nrow = 16, ncol = 24)
  
  # Fill the middle 22 columns with the matrix data
  mat[, 2:23] <- mat_data
  
  # Set row and column names
  rownames(mat) <- LETTERS[1:16]
  colnames(mat) <- 1:24
  
  return(mat)
}

plate_matrix <- ConvertMatrixtoPlate(RperU_values)

svg("ReadsperUMI.svg", width = 40,  height = 10) 
barplot(RperU_df$RperU_matrix, 
        names.arg = RperU_df$Sample, 
        col = "pink", 
        main = "Barplot of RperU_matrix values for each Sample", 
        xlab = "Sample", 
        ylab = "ReadsPerUMI", 
        las = 2)
        
dev.off()
totalreads <- reads[1,]
totalreads <- as.data.frame(totalreads)
colnames(totalreads) <- colnames(reads) 
totalreads_matrix <- as.matrix(totalreads)
head(totalreads_matrix) 
totalreads_matrix <- totalreads_matrix[, match(metadata$well, colnames(totalreads_matrix))] 
head(totalreads_matrix) 
totalreads_df <- as.data.frame(totalreads_matrix) 
head(totalreads_df) 
totalreads_df$Sample <- rownames(totalreads_df)


head(totalreads_df)


genes_detect <- colSums(A2>0)
genes_detect <- as.data.frame(genes_detect)
genes_detect$Sample <- rownames(genes_detect) 
genes_detect <- genes_detect[match(metadata$well, genes_detect$Sample), ]
UMI <- colSums(A2)
UMI <- as.data.frame(UMI)
UMI$Sample <- rownames(UMI)
UMI <- UMI[match(metadata$well, UMI$Sample), ]

str(UMI)

str(genes_detect) 

str(totalreads_df)

totalreads_df$reads_per_gene <- totalreads_df$totalreads_matrix / genes_detect$genes_detect
totalreads_df$reads_per_UMI <- totalreads_df$totalreads_matrix / UMI$UMI
totalreads_df$UMI_per_gene <- UMI$UMI / genes_detect$genes_detect

## Plotting 
svg ("Reads_per_genes_detected.svg", width = 15, height = 25)
plot(totalreads_df$totalreads_matrix, genes_detect$genes_detect, 
     xlab = "Total Reads", 
     ylab = "Number of Genes Detected", 
     main = "Reads Per Genes Detected", 
     pch = 19,   # Points as filled circles
     col = "blue", 
     cex = 0.7)
dev.off()



svg ("Reads_per_UMI_detected.svg", width = 15, height = 25)
plot(totalreads_df$totalreads_matrix, UMI$UMI, 
     xlab = "Total Reads", 
     ylab = "Number of UMIs detected",
     main = "Reads Per UMI detected", 
     pch = 19,   # Points as filled circles
     col = "hotpink", 
     cex = 0.7)
dev.off()

write.csv(totalreads_df, "totalreads_df.csv")

# Check the result
head(totalreads_df)
## UMI Counts Log Barplot Across Samples 
svg("UMICounts_log_barplot.svg", width = 10, height = 7.5)  
par(mar = c(20, 5, 5, 2))  


barplot(log(colSums(A2)),
        main = "UMI Counts Across Samples", 
        ylab = "log(colSums(A))",
        las = 2,  
        cex.names = 1.0,  
        col = "lightpink")  
dev.off()

## Genes Across Samples 

svg("Genes_log_barplot.svg", width = 10, height = 7.5)  
par(mar = c(20, 5, 5, 2))  


barplot(log(colSums(A>0)), 
        main = "Genes Across Samples", 
        ylab = "log(colSums(A>0))",
        las = 2,  
        cex.names = 1.0,  
        col = "orange")  
dev.off()

pseudocount = 100
B2 = A[,colSums(A) >= 5000]
#Cutoff of 10,000 UMIs has pass of 81% (26/32) 
#Cutoff of 1000 gives all 32
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,]

library(ComplexHeatmap)

# drop, 14s, 15s, 18s, 32s 
cor_matrix = cor(B4, method = "pearson")
Heatmap(cor_matrix)
svg("CorrelationHeatmap_TS_RC1November2024.svg", width = 24, height = 24)
Heatmap(cor_matrix)
dev.off()