### 1/6/2025 

#### August Data, R and C1 


# Instead of combining experimental treatments keep them separate before running DESeq2 

setwd("C:/Users/taylo/Desktop/AugustSequencingNelms2024/")
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

colnames(A) <- sub(".*_(\\d+s)\\.bam\\.tsv$", "\\1", colnames(A)) # Renaming columns to just make them the number and s 
colnames(A) # just check what everything is 
A[is.na(A)] = 0
A = A[rowSums(A) > 0,]
rownames(A) = annots[rownames(A)]
dim(A)
[1] 34747    96
A <- A[grepl("^Zm", rownames(A)), ]
dim(A)
[1] 26257    96

## check for duplicated rownames with duplicated(rownames(A)) in base R studio 
# can also use A[duplicated(rownames(A))] to see the exact positons of duplicates 

for (g in unique(rownames(A)[duplicated(rownames(A))])) {
    i = which(rownames(A) == g)
    A[i[1],] = colSums(A[i,])
    A = A[-i[-1],]
}


A = A[,-10] #remove 19s as it was strange 
A = A[,order(as.numeric(sub('s','',colnames(A))))]
colnames(A) 
A = A[,-c(32:96)] #removing columns 32 through 96 
colnames(A)

samples <- c("5s", "6s", "7s", "8s", "25s", "26s", "27s", "28s",   # R
             "17s", "18s", "20s", "21s", "22s", "23s", "24s",       # C1
             "1s", "2s", "3s", "4s", "9s", "10s", "11s", "12s",     # R+C1
             "13s", "14s", "15s", "16s", "29s", "30s", "31s", "32s") # Control
conditions <- c(rep("R", 8),       # 8 samples for R
                rep("C1", 7),     # 7 samples for C1
                rep("R+C1", 8),   # 8 samples for R+C1
                rep("Control", 8)) # 8 samples for Control
colData <- data.frame(
  condition = factor(conditions, levels = c("R", "C1", "R+C1", "Control"))
)

rownames(colData) <- samples            
counts <- A[,samples] ## now re-ordered based on the samples that we gave above 
counts[is.na(counts)] <- 0 #replaced NA with 0 

pseudocount = 100
B2 = A[,colSums(A) >= 5000]
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,] ### Values in B4 are now TPM normalized with pseudocount of 100 and log transformed
logTPM <- t(scale(t(B4)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)

R_control <- results(dds, contrast = c("condition", "R", "Control"))
C1_control <- results(dds, contrast = c("condition", "C1", "Control"))
R_C1_control <- results(dds, contrast = c("condition", "R+C1", "Control"))
summary (R_control) 
summary(C1_control)
summary(R_C1_control)

## summary (R_control) 

out of 26098 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4948, 19%
LFC < 0 (down)     : 3826, 15%
outliers [1]       : 0, 0%
low counts [2]     : 8099, 31%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

## summary(C1_control)

out of 26098 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1245, 4.8%
LFC < 0 (down)     : 1822, 7%
outliers [1]       : 0, 0%
low counts [2]     : 11134, 43%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

## summary(R_C1_control)

out of 26098 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4207, 16%
LFC < 0 (down)     : 3212, 12%
outliers [1]       : 0, 0%
low counts [2]     : 8605, 33%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results



### Important Gene Names 
GRMZM2G018724	C2
GRMZM2G701063 pl1 
GRMZM2G016241 bz2 
GRMZM2G005066 C1
GRMZM2G026930 A1
GRMZM2G165390 bz1 
GRMZM2G335358 p1
Zm00001eb429330  R1 

C2 Zm00001eb373660 
pl1 Zm00001eb048110
bz2 Zm00001eb159020
c1 Zm00001eb374220	Zm00001eb374230
A1 Zm00001eb014420	Zm00001eb014430
p1 Zm00001eb278680
R1 Zm00001eb429330  
--- C2, pl1, bz2, c1, A1, bz1, p1, and R1 

## KEGG Analysis -- R_C1_Control 


## Anthocyanin Genes Enriched in RC1 Sample 

R_C1_control_filtered <- R_C1_control[!is.na(R_C1_control$log2FoldChange) & !is.na(R_C1_control$padj), ]
Res_RC1 <- R_C1_control_filtered[(R_C1_control_filtered$log2FoldChange > 1) & (R_C1_control_filtered$padj < 0.001), ]

Res_RC1$FCRank = rank(-Res_RC1[,2])
Res_RC1$pRank = rank(Res_RC1$pvalue)

#Res_RC1_sorted <- Res_RC1[order(Res_RC1$padj), ]
#RC1_top_100_significant_genes <- Res_RC1_sorted[1:100, ]
#head(RC1_top_100_significant_genes)
#write.csv(RC1_top_100_significant_genes, file = "top_100_significant_genes_RC1.csv", row.names = TRUE)

anthGenes = c('Zm00001eb198030', 'Zm00001eb062510', 'Zm00001eb067380', 'Zm00001eb159020', 'Zm00001eb229190', 'Zm00001eb374230', 'Zm00001eb048110')
names(anthGenes) = c('C2', 'CHI1', 'F3H', 'A1', 'A2','bz1', 'bz2')

Anthocyanin_RC1 <- Res_RC1[rownames(Res_RC1) %in% anthGenes, ]
rownames(Anthocyanin_RC1) <- names(anthGenes)[match(rownames(Anthocyanin_RC1), anthGenes)]

Wald test p-value: condition R+C1 vs Control 
DataFrame with 7 rows and 8 columns
      baseMean log2FoldChange     lfcSE      stat       pvalue         padj
     <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
bz2    59.7875        5.83583  0.287055  20.32999  6.98157e-92  8.72646e-89
CHI1   66.6688        6.03126  0.293442  20.55349  7.16261e-94  9.64143e-91
F3H     4.6588        5.64611  0.746705   7.56137  3.98852e-14  7.17319e-13
A1     81.2903        5.60831  0.351669  15.94772  2.95485e-57  9.23339e-55
C2    545.1632        6.11186  0.260135  23.49497 4.59176e-122 1.60702e-118
A2    159.4353        5.91214  0.383079  15.43323  9.78431e-54  2.34542e-51
bz1    13.7483        6.02132  0.561248  10.72845  7.48409e-27  3.79606e-25
        FCRank     pRank
     <numeric> <numeric>
bz2         57        13
CHI1        44        12
F3H         73       428
A1          76        26
C2          39         5
A2          51        31
bz1         45       125

## KEGG RC1 
library(clusterProfiler)
library(pathview)
search_kegg_organism('zma', by = 'kegg_code')
    kegg_code scientific_name common_name
678       zma        Zea mays       maize

RC1 <- read.csv("R_C1_ENTREZ.csv")
RC1_cleaned <- RC1[RC1$GenBank.Gene != "", ]
RC1_cleaned <- na.omit(RC1_cleaned)
RC1_gene <- RC1_cleaned$GenBank.Gene
RC1_gene <- gsub("^LOC", "", RC1_gene)


kegg_results_RC1 <- enrichKEGG(gene = RC1_gene,
                           organism = "zma",    # KEGG code for maize
                           pvalueCutoff = 0.05)

head(kegg_results_RC1)

dotplot(kegg_results_RC1, showCategory = 5)






## Anthocyanin Genes Enriched in R_Control Sample 
R_control_filtered <- R_control[!is.na(R_control$log2FoldChange) & !is.na(R_control$padj), ]
Res_R <- R_control_filtered[(R_control_filtered$log2FoldChange > 1) & (R_control_filtered$padj < 0.001), ]

Res_R$FCRank = rank(-Res_R[,2])
Res_R$pRank = rank(Res_R$pvalue)


anthGenes = c('Zm00001eb198030', 'Zm00001eb062510', 'Zm00001eb067380', 'Zm00001eb159020', 'Zm00001eb229190', 'Zm00001eb374230', 'Zm00001eb048110')
names(anthGenes) = c('C2', 'CHI1', 'F3H', 'A1', 'A2','bz1', 'bz2')

Anthocyanin_R <- Res_R[rownames(Res_R) %in% anthGenes, ]
rownames(Anthocyanin_R) <- names(anthGenes)[match(rownames(Anthocyanin_R), anthGenes)]
head(Anthocyanin_R)

Wald test p-value: condition R vs Control 
DataFrame with 2 rows and 8 columns
     baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
bz2   59.7875        1.93900  0.306400   6.32832 2.47848e-10 2.53984e-09
bz1   13.7483        2.26017  0.598507   3.77634 1.59149e-04 6.58426e-04
       FCRank     pRank
    <numeric> <numeric>
bz2       798       661
bz1       578      1803

## Anthocyanin Genes Enriched in C1_Control Sample 

C1_control_filtered <- C1_control[!is.na(C1_control$log2FoldChange) & !is.na(C1_control$padj), ]
Res_C1 <- C1_control_filtered[(C1_control_filtered$log2FoldChange > 1) & (C1_control_filtered$padj < 0.001), ]

Res_C1$FCRank = rank(-Res_C1[,2])
Res_C1$pRank = rank(Res_C1$pvalue)

anthGenes = c('Zm00001eb198030', 'Zm00001eb062510', 'Zm00001eb067380', 'Zm00001eb159020', 'Zm00001eb229190', 'Zm00001eb374230', 'Zm00001eb048110')
names(anthGenes) = c('C2', 'CHI1', 'F3H', 'A1', 'A2','bz1', 'bz2')

Anthocyanin_C1 <- Res_C1[rownames(Res_C1)%in% anthGenes, ]
rownames(Anthocyanin_C1) <- names(anthGenes)[match(rownames(Anthocyanin_C1), anthGenes)]
head(Anthocyanin_C1)

log2 fold change (MLE): condition C1 vs Control 
Wald test p-value: condition C1 vs Control
DataFrame with 4 rows and 8 columns
     baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
bz2   59.7875        2.54358  0.305509   8.32573 8.38106e-17 1.60852e-14
C2   545.1632        1.94762  0.272502   7.14717 8.85843e-13 8.30525e-11
A2   159.4353        1.67462  0.404859   4.13631 3.52932e-05 5.39672e-04
bz1   13.7483        2.77222  0.592453   4.67922 2.87974e-06 6.08035e-05
       FCRank     pRank
    <numeric> <numeric>
bz2        14        18
C2         30        32
A2         50       168
bz1         8       121