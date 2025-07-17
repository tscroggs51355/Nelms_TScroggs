# RC1 Novemeber 2024 QC analysis 

setwd("C:/Users/taylo/Desktop/RC1_November2024/")
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

colnames(A) <- gsub(".*_(\\d+s)\\.bam\\.tsv", "\\1", colnames(A))
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

R Alone: 1s - 4s 
R + Reporter: 5s - 8s 
R + C1: 9s - 12s 
C1 Alone: 13s - 16s 
R + C1 + Reporter: 17s - 20s 
C1 + Reporter: 21s - 24s 
Reporter: 25s - 28s 
mCherry Control: 29s - 32s 

samples <- c("1s", "2s", "3s", "4s", "5s", "6s", "7s", "8s",   # R1 Alone
             "13s", "14s", "15s", "16s", "21s", "22s", "23s", '24s',      # C1 Alone
             "9s", "10s", "11s", "12s", "17s", "18s", "19s", "20s",     # R+C1
             "25s", "26s", "27s", "28s", "29s", "30s", "31s", "32s") # Control
conditions <- c(rep("R1", 8),       # 8 samples for R
                rep("C1", 8),     # 7 samples for C1
                rep("R1+C1", 8),   # 8 samples for R+C1
                rep("Control", 8)) # 8 samples for Control
names(conditions) = samples
colData <- data.frame(
  condition = factor(conditions, levels = c("R1", "C1", "R1+C1", "Control"))
)

counts <- A[,samples] ## now re-ordered based on the samples that we gave above 
counts[is.na(counts)] <- 0 #replaced NA with 0 

pseudocount = 100
B2 = A[,colSums(A) >= 5000]
TPM = sweep(B2,2,colSums(B2),'/')*10^6
B3 = log(TPM + pseudocount, 10)
logTPM <- t(scale(t(B3)))


library(DESeq2)
rownames(colData) <- samples            

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

R1_control <- results(dds, contrast = c("condition", "R1", "Control"))
C1_control <- results(dds, contrast = c("condition", "C1", "Control"))
R1_C1_control <- results(dds, contrast = c("condition", "R1+C1", "Control"))
summary (R1_control) 
summary(C1_control)
summary(R1_C1_control)

PVal_R1 <- R1_control$pvalue
PVal_C1 <- C1_control$pvalue
PVal_R1_C1 <- R1_C1_control$pvalue

combined_pvalues <- data.frame(
  R1= PVal_R1,
  C1= PVal_C1,
  R1_C1= PVal_R1_C1
)
rownames(combined_pvalues) = rownames(R1_control)
X2 = -2*rowSums(log(combined_pvalues))
combined_pvalues$p_pool = pchisq(X2, 2*ncol(combined_pvalues), lower.tail=F) + dchisq(X2,2*3)
Pooled_pvalue <- combined_pvalues$p_pool
Adjusted_pvalue <- p.adjust(Pooled_pvalue, method = "holm") 
significant_adjustedpvalue <- sum(Adjusted_pvalue <= 0.05, na.rm = TRUE)
significant_adjustedpvalue

combined_pvalues$padj = p.adjust(combined_pvalues$p_pool, method = "holm") 

combined_pvalues$maxLog2 = apply(cbind(R1_control[,2], C1_control[,2], R1_C1_control[,2]), 1, function(xx) { xx[rank(-abs(xx), ties.method = 'first') == 1] })

combined_pvalues$Significant = (combined_pvalues$padj <= .05) & (abs(combined_pvalues$maxLog2) >= 1)
combined_pvalues$Significant[is.na(combined_pvalues$Significant)] = FALSE

combined_pvalues = combined_pvalues[order(-combined_pvalues$Significant, combined_pvalues$padj),]
sigGenes = rownames(combined_pvalues)[combined_pvalues$Significant]

## Figure 3

set.seed(1)
library(ComplexHeatmap)

valid_samples <- intersect(samples, colnames(logTPM))
hmMat <- logTPM[
  which(rank(R1_C1_control$pvalue) <= 100),
  valid_samples
]

names(conditions) <- valid_samples

valid_group_labels <- conditions[valid_samples]


column_anno <- HeatmapAnnotation(
  Group = valid_group_labels,
  col = list(Group = c(
    "R1"      = "#FF9999",
    "C1"      = "#99CCFF",
    "R1+C1"    = "#66FF66",
    "Control" = "#CCCCCC"
  )),
  show_legend = TRUE
)


hm <- Heatmap(
  hmMat,
  top_annotation = column_anno,
  show_row_names = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_distance_rows = "euclidean",
  cluster_columns = FALSE,
  show_column_names = TRUE
)

svg("Fig3_HeatmapwithAnnotation_November2024Data.svg", width = 10, height = 20) 
plot(hm)
dev.off()










#########################################################
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1561   22184   32576   46408   64799  144264


summary(colSums(A>0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    969    6974    7958    7682    9799   12359
dim(A)
[1] 21024    32

 colnames(A)
 [1] "10s" "11s" "12s" "13s" "14s" "15s" "16s" "17s" "18s" "19s" "1s"  "20s"
[13] "21s" "22s" "23s" "24s" "25s" "26s" "27s" "28s" "29s" "2s"  "30s" "31s"
[25] "32s" "3s"  "4s"  "5s"  "6s"  "7s"  "8s"  "9s"

R Alone: 1s - 4s 
R + Reporter: 5s - 8s 
R + C1: 9s - 12s 
C1 Alone: 13s - 16s 
R + C1 + Reporter: 17s - 20s 
C1 + Reporter: 21s - 24s 
Reporter: 25s - 28s 
mCherry Control: 29s - 32s 
####

Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/RC1_November2024/Mapped_Data/stringtie_out/read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)

colnames(reads) <- sub(".*_(\\d+s)\\.bam$", "\\1", colnames(reads))

## Getting UMI Counts 
RperU = (reads[1,])/(colSums(A))
write.csv(RperU, "ReadsperUMICalculation_RC1November2024Sequencing.csv")

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

