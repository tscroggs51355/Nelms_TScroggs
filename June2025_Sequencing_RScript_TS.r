## Moving to R Studio Now 

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



#########################################################
> summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   2478   65092  120211  163502  251935  49473

>  summary(colSums(A > 0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1132    9220   10903   10972   13886   15910
####

Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/June2025_Sequencing/Mapped_Data/stringtie_out/read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- reads[, !grepl("_unsorted\\.bam$", colnames(reads))]
colnames(reads) <- sub("^Mapped_Data\\.hisat2_out\\.Justin_", "", colnames(reads))
colnames(reads) <- sub("\\.bam$", "", colnames(reads))

## Getting UMI Counts 
RperU = (reads[1,])/(colSums(A))
write.csv(RperU, "ReadsperUMICalculation_June2025Sequencing.csv")

## UMI Counts Log Barplot Across Samples 
svg("UMICounts_log_barplot.svg", width = 10, height = 7.5)  
par(mar = c(20, 5, 5, 2))  


barplot(log(colSums(A)), 
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
B2 = A[,colSums(A) >= 2000]
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,]

library(ComplexHeatmap)

cor_matrix_2 <- B4[, !colnames(B4) %in% c("STP2_S3_85s", "STP1_S1_85s")]
cor_matrix = cor(B4, method = "pearson")
cor_matrix_2 = cor(cor_matrix_2, method = "pearson")
Heatmap(cor_matrix)
svg("CorrelationHeatmap_TS_June2025Sequencing.svg", width = 24, height = 24)
Heatmap(cor_matrix)
dev.off()
svg("Correlation_2.svg", width = 24, height = 24)
Heatmap(cor_matrix_2)
dev.off()


pairs(unique(A[,c('STP1_S1_89s','STP1_S1_90s','STP2_S3_85s','STP2_S3_86s')]) + 1, log = 'xy', pch = 19, cex = .8)

##No DNaseI, No Delay 
#svg("NoDNaseI_NoDelay.svg", width = 24, height =24)
png("NoDNaseI_NoDelay.png", width = 24, height =24)
pairs(unique(A[,c('STP1_S1_81s', 'STP1_S1_82s', 'STP1_S1_83s', 'STP1_S1_84s', 'STP2_S3_90s', 'STP2_S3_91s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X1 <- A[,c('STP1_S1_81s', 'STP1_S1_82s', 'STP1_S1_83s', 'STP1_S1_84s', 'STP2_S3_90s', 'STP2_S3_91s')]

#STP1_S1_81s
#STP1_S1_82s
#STP1_S1_83s 
#STP1_S1_84s
#STP2_S3_90s
#STP2_S3_91s

#No DNase I, 1 hr Delay 
#svg("NoDNaseI_1hrDelay.svg", width = 24, height =24)
png("NoDNaseI_1hrDelay.png", width = 24, height =24)
pairs(unique(A[,c('STP1_S1_85s', 'STP1_S1_86s', 'STP2_S3_92s', 'STP2_S3_93s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X2 <- A[,c('STP1_S1_85s', 'STP1_S1_86s', 'STP2_S3_92s', 'STP2_S3_93s')]

#STP1_S1_85s
#STP1_S1_86s
#STP2_S3_92s
#STP2_S3_93s

#No DnaseI, 4 hr delay 
#svg("NoDNaseI_4hrDelay.svg", width = 24, height =24)
png("NoDNaseI_4hrDelay.png", width = 24, height =24)
pairs(unique(A[,c('STP2_S3_94s', 'STP2_S3_95s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X3 <- A[,c('STP2_S3_94s', 'STP2_S3_95s')]

#STP2_S3_94s
#STP2_S3_95s

#Dnase I no Delay
#svg("DNaseI_NoDelay.svg", width = 24, height = 24) 
png("DNaseI_NoDelay.png", width = 24, height = 24) 
pairs(unique(A[,c('STP1_S1_87s', 'STP1_S1_88s','STP2_S3_81s', 'STP2_S3_82s', 'STP2_S3_83s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X4 <- A[,c('STP1_S1_87s', 'STP1_S1_88s','STP2_S3_81s', 'STP2_S3_82s', 'STP2_S3_83s')]

#STP1_S1_87s
#STP1_S1_88s
#STP2_S3_81s
#STP2_S3_82s
#STP2_S3_83s

#DNase I 1 hr delay 
#svg("DNaseI_1hrDelay.svg", width = 24, height = 24) 
png("DNaseI_1hrDelay.png", width = 24, height = 24) 
pairs(unique(A[,c('STP1_S1_89s', 'STP1_S1_90s','STP2_S3_84s', 'STP2_S3_85s', 'STP2_S3_86s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X5 <- A[,c('STP1_S1_89s', 'STP1_S1_90s','STP2_S3_84s', 'STP2_S3_85s', 'STP2_S3_86s')]
#STP1_S1_89s
#STP1_S1_90s
#STP2_S3_84s
#STP2_S3_85s
#STP2_S3_86s

#DNaseI, 1 hr delay, Double Volume, and Source well(STP1_S1_90s)
#svg("DNaseI_1hrdelay_double_volumeRNA_source.svg", width = 24, height = 24) 
png("DNaseI_1hrdelay_double_volumeRNA_source.png", width = 24, height = 24) 
pairs(unique(A[,c('STP1_S1_1s', 'DP1_10_S2_1s', 'STP1_S1_90s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
X6 <- A[,c('STP1_S1_1s', 'DP1_10_S2_1s', 'STP1_S1_90s')]
#STP1_S1_1s
#DP1_10_1s 
#STP1_S1_90s

#DNaseI, 1 hr delay, Quadruplicate Volume, and Source Well(STP2_S3_85s)
#svg("DNaseI_1hrdelay_quad_volumeRNA_source.svg", width = 24, height = 24) 
png("DNaseI_1hrdelay_quad_volumeRNA_source.png", width = 24, height = 24) 
pairs(unique(A[,c('STP2_S3_1s', 'FP2_5_S4_1s', 'STP2_S3_85s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
A7 <- A[,c('STP2_S3_1s', 'FP2_5_S4_1s', 'STP2_S3_85s')]
#STP2_S3_1s
#FP2_5_1s 
#STP2_S3_85s

#DNaseI, 4 hr delay, 
#svg("DNaseI_4hrdelay.svg", width = 24, height = 24) 
png("DNaseI_4hrdelay.png", width = 24, height = 24) 
pairs(unique(A[,c('STP2_S3_87s', 'STP2_S3_88s', 'STP2_S3_89s')]) + 1, log = 'xy', pch = 19, cex = .8)
dev.off()
A8 <- A[,c('STP2_S3_87s', 'STP2_S3_88s', 'STP2_S3_89s')]
#STP2_S3_87s
#STP2_S3_88s
#STP2_S3_89s

##No DNaseI, No Delay 
X1 <- A[,c('STP1_S1_81s', 'STP1_S1_82s', 'STP1_S1_83s', 'STP1_S1_84s', 'STP2_S3_90s', 'STP2_S3_91s')]
#No DNase I, 1 hr Delay 
X2 <- A[,c('STP1_S1_85s', 'STP1_S1_86s', 'STP2_S3_92s', 'STP2_S3_93s')]
#No DnaseI, 4 hr delay 
X3 <- A[,c('STP2_S3_94s', 'STP2_S3_95s')]
#Dnase I no Delay
X4 <- A[,c('STP1_S1_87s', 'STP1_S1_88s','STP2_S3_81s', 'STP2_S3_82s', 'STP2_S3_83s')]
#DNase I 1 hr delay 
X5 <- A[,c('STP1_S1_89s', 'STP1_S1_90s','STP2_S3_84s', 'STP2_S3_85s', 'STP2_S3_86s')]
#DNaseI, 1 hr delay, Double Volume, and Source well(STP1_S1_90s)
X6 <- A[,c('STP1_S1_1s', 'DP1_10_S2_1s', 'STP1_S1_90s')]
#DNaseI, 1 hr delay, Quadruplicate Volume, and Source Well(STP2_S3_85s)
X7 <- A[,c('STP2_S3_1s', 'FP2_5_S4_1s', 'STP2_S3_85s')]
#DNaseI, 4 hr delay, 
X8 <- A[,c('STP2_S3_87s', 'STP2_S3_88s', 'STP2_S3_89s')]

pseudocount = 100
B2 = A[,colSums(A) >= 2000]
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,]
data_list <- list(X1, X2, X3, X4, X5, X6, X7, X8)
col_sums <- sapply(data_list, colSums)
col_sums_genes <- sapply(data_list, function(x) colSums(x) > 0)
boxplot(col_sums, ylim = c(0000,500000), main = "UMI Counts Across Samples", col = "pink")
boxplot(col_sums_genes, ylim = c(0000,20000), main = "Genes Per Sample", col = "lightgreen")


pseudocount <- 100

processed_list <- lapply(list(X1, X2, X3, X4, X5, X6, X7, X8), function(A) {
  B2 <- A[, colSums(A) >= 2000]
  B3 <- log(sweep(B2, 2, colSums(B2), '/') * 10^6 + pseudocount, 10)
  B4 <- B3[rowSums(B2 >= 10) >= 2, ]
  return(B4)
})

B4_X1 <- processed_list[[1]]
B4_X2 <- processed_list[[2]]
B4_X3 <- processed_list[[3]]
B4_X4 <- processed_list[[4]]
B4_X5 <- processed_list[[5]]
B4_X6 <- processed_list[[6]]
B4_X7 <- processed_list[[7]]
B4_X8 <- processed_list[[8]]

svg("X1_correlation.svg", width = 12, height = 8)
cor_matrix_X1 = cor(B4_X1, method = "pearson")
Heatmap(cor_matrix_X1)
dev.off()

svg("X2_correlation.svg", width = 12, height = 8)
cor_matrix_X2 = cor(B4_X2, method = "pearson")
Heatmap(cor_matrix_2)
dev.off()

svg("X3_correlation.svg", width = 12, height = 8)
cor_matrix_X3 = cor(B4_X3, method = "pearson")
Heatmap(cor_matrix_X3)
dev.off()

svg("X4_correlation.svg", width = 12, height = 8)
cor_matrix_X4 = cor(B4_X4, method = "pearson")
Heatmap(cor_matrix_X4)
dev.off()

svg("X5_correlation.svg", width = 12, height = 8)
cor_matrix_X5 = cor(B4_X5, method = "pearson")
Heatmap(cor_matrix_X5)
dev.off()

svg("X6_correlation.svg", width = 12, height = 8)
cor_matrix_X6 = cor(B4_X6, method = "pearson")
Heatmap(cor_matrix_X6)
dev.off()

svg("X7_correlation.svg", width = 12, height = 8)
cor_matrix_X7 = cor(B4_X7, method = "pearson")
Heatmap(cor_matrix_X7)
dev.off()

svg("X8_correlation.svg", width = 12, height = 8)
cor_matrix_X8 = cor(B4_X8, method = "pearson")
Heatmap(cor_matrix_X8)
dev.off()

STP1 <- A[,c('DP1_10_S2_1s', 'STP1_S1_83s', 'STP1_S1_84s', 'STP1_S1_85s', 'STP1_S1_86s', 'STP1_S1_87s', 'STP1_S1_88s', 'STP1_S1_89s', 'STP1_S1_90s', 'STP1_S1_1s', 'STP1_S1_81s', 'STP1_S1_82s')]
write.csv(as.data.frame(STP1), "STP1_UMICounts.csv", row.names = TRUE)

STP2 <- A[, c('FP2_5_S4_1s', 'STP2_S3_1s', 'STP2_S3_81s', 'STP2_S3_82s', 'STP2_S3_83s','STP2_S3_84s', 'STP2_S3_85s', 'STP2_S3_86s', 'STP2_S3_87s', 'STP2_S3_88s', 'STP2_S3_89s', 'STP2_S3_90s', 'STP2_S3_91s', 'STP2_S3_92s', 'STP2_S3_93s', 'STP2_S3_94s', 'STP2_S3_95s')]
write.csv(as.data.frame(STP2), "STP2_UMICounts.csv", row.names = TRUE)

data_list_1 <- list(STP1, STP2)
col_sums_Batch <- sapply(data_list_1, colSums)
col_sums_genes_Batch <- sapply(data_list_1, function(x) colSums(x > 0))
boxplot(col_sums_Batch, ylim = c(0000,500000), main = "UMI Counts Across Samples Per Batch", col = "lightpink")
write.csv(as.data.frame(col_sums_Batch), "col_sums_Batch.csv", row.names = TRUE)

boxplot(col_sums_genes_Batch, ylim = c(0000,20000), main = "Genes Per Sample Per Batchy", col = "lightgreen")
write.csv(as.data.frame(col_sums_genes_Batch), "col_sums_genes_Batch.csv", row.names = TRUE)


col_sums_Batch_df <- do.call(rbind, lapply(col_sums_Batch, as.data.frame))
write.csv(col_sums_Batch_df, "col_sums_Batch.csv", row.names = TRUE)

col_sums_genes_Batch_df <- do.call(rbind, lapply(col_sums_genes_Batch, as.data.frame))
write.csv(col_sums_genes_Batch_df, "col_sums_genes_Batch.csv", row.names = TRUE)

