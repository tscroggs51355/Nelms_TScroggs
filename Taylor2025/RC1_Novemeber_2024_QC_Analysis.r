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

#ColSums Less than 5000, true for 14s, 15s, 18s, 32s 

samples <- c("1s", "2s", "3s", "4s", "5s", "6s", "7s", "8s",   # R1 Alone
             "13s",  "16s", "21s", "22s", "23s", '24s',      # C1 Alone
             "9s", "10s", "11s", "12s", "17s", "19s", "20s",     # R1+C1
             "25s", "26s", "27s", "28s", "29s", "30s", "31s") # Control
conditions <- c(rep("R1", 8),       # 8 samples for R1
                rep("C1", 6),     # 7 samples for C1
                rep("R1+C1", 7),   # 8 samples for R+C1
                rep("Control", 7)) # 8 samples for Control
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
#547
combined_pvalues$padj = p.adjust(combined_pvalues$p_pool, method = "holm") 

combined_pvalues$maxLog2 = apply(cbind(R1_control[,2], C1_control[,2], R1_C1_control[,2]), 1, function(xx) { xx[rank(-abs(xx), ties.method = 'first') == 1] })

combined_pvalues$Significant = (combined_pvalues$padj <= .05) & (abs(combined_pvalues$maxLog2) >= 1)
combined_pvalues$Significant[is.na(combined_pvalues$Significant)] = FALSE

combined_pvalues = combined_pvalues[order(-combined_pvalues$Significant, combined_pvalues$padj),]
sigGenes = rownames(combined_pvalues)[combined_pvalues$Significant]

## Figure 3

set.seed(1)
library(ComplexHeatmap)
library(grid)
anthGenes = c('C2', 'CHI1', 'CHI3', 'CHI4', 'CHI5', 'CHI6', 'F3H', 'A1', 'A2', 'bz1', 'bz2')
names(anthGenes) = c('Zm00001eb198030', 'Zm00001eb062510', 'Zm00001eb211840',
 'Zm00001eb256920', 'Zm00001eb238460', 'Zm00001eb163860', 'Zm00001eb067380'#F3H, 
 'Zm00001eb159020', 'Zm00001eb229190', 'Zm00001eb374230', 'Zm00001eb048110')



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
    "R1"      = "#F1BB7B",
    "C1"      = "#FD6467",
    "R1+C1"    = "#5B1A18",
    "Control" = "#D6C6B9"
  )),
  show_legend = TRUE
)

heatmap_rows <- rownames(hmMat)
label_colors <- rep("black", length(heatmap_rows))
names(label_colors) <- heatmap_rows
label_colors[names(anthGenes)[names(anthGenes) %in% heatmap_rows]] <- "red"

# Plot heatmap with colored row labels
hm <- Heatmap(
  hmMat,
  top_annotation = column_anno,
  show_row_names = TRUE,
  row_names_gp = gpar(col = label_colors),
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

hmGenes <- rownames(hmMat)[row_order(hm)]

anthGenes = c('C2', 'CHI1', 'CHI3', 'CHI4', 'CHI5', 'CHI6', 'F3H', 'A1', 'A2', 'bz1', 'bz2')
names(anthGenes) = c('Zm00001eb198030', 'Zm00001eb062510', 'Zm00001eb211840', 'Zm00001eb256920', 'Zm00001eb238460', 'Zm00001eb163860', 'Zm00001eb067380', 'Zm00001eb159020', 'Zm00001eb229190', 'Zm00001eb374230', 'Zm00001eb048110')
anthGenes = anthGenes[-c(3:5)]

cbind(1:100, anthGenes[hmGenes])

library(reshape2)
library(ggplot2)
anthExp = TPM[names(anthGenes),]
rownames(anthExp) = anthGenes
colnames(anthExp) = conditions[colnames(anthExp)]
Antho_Melt <- melt(anthExp)

svg("Fig3_AnthGenesBarplot.svg", width = 9, height = 18) 

library(ggplot2)
library(wesanderson)

Fig3_AnthGenes <- ggplot(Antho_Melt, aes(x = Var2, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~Var1, scales = "free_y", nrow = 8)

plot(Fig3_AnthGenes)

dev.off()

## 




###--- Load and process GAF file ---###
gaf <- read.delim("1.1_GOMAP-output.gaf.gz", header = FALSE, comment.char = "!")
gaf$Gene <- sub("_T[0-9]+$", "", gaf$V2)  # Collapse isoform IDs to gene-level
gene2go <- split(gaf$V5, gaf$Gene)
gene2go <- lapply(gene2go, unique)

###--- Prepare gene list for topGO ---###
geneUniverse <- names(sigGenes)
geneList <- factor(as.integer(geneUniverse %in% sigGenes))
names(geneList) <- geneUniverse

###--- Run topGO enrichment analysis ---###
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = function(p) p == 1,
              annot = annFUN.gene2GO,
              gene2GO = gene2go)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
goResults <- GenTable(GOdata, classicFisher = resultFisher,
                      topNodes = length(score(resultFisher)))
write.csv(goResults, "GO_enrichment_topGO_B73v5.csv", row.names = FALSE)



### 7/30/2025
#--- Load libraries ---#
library(topGO)
library(dplyr)

#--- Load and process GAF file ---#
gaf <- read.delim("1.1_GOMAP-output.gaf.gz", header = FALSE, comment.char = "!")
gaf$Gene <- sub("_T[0-9]+$", "", gaf$V2)  # Collapse isoform IDs to gene-level
gene2go <- split(gaf$V5, gaf$Gene)
gene2go <- lapply(gene2go, unique)

#--- Define universe ---#
geneUniverse <- names(gene2go)



##Clustering, Kmeans, SigGenes 
## Reminder - sigGenes = rownames(combined_pvalues)[combined_pvalues$Significant]

# Apply k-means clustering with k = 5
set.seed(1)
km_res <- kmeans(logTPM[sigGenes, ], centers = 5)

# Create cluster assignment data.frame
cluster_df <- data.frame(Gene = sigGenes, Cluster = km_res$cluster)



## Figuring out why kmeans isn't running -- 

> sum(is.na(logTPM))       # Count of NA
[1] 2464
> sum(is.nan(logTPM))      # Count of NaN
[1] 2464
> sum(is.infinite(logTPM)) # Count of infinite 
[1] 0

# bc there was this error 

> km_res <- kmeans(logTPM, centers = 5)
Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)

##Cleaned up logTPM to remove NA values 

logTPM<- logTPM[apply(logTPM, 1, function(x) all(is.finite(x))), ]







#--- Define cluster assignments ---#
cluster_df <- data.frame(Gene = sigGenes, Cluster = clusters)


set.seed(1)
library(ComplexHeatmap)

hm <- Heatmap(logTPM[sigGenes[1:1000],], 
        show_row_names = FALSE,  
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "euclidean",
        cluster_rows = TRUE, 
        row_split = 5, 
        show_column_names = TRUE)

plot(hm)

clusterlist <- row_order(hm)
Clusters <- lapply(clusterlist, function(xx) { rownames(combined_pvalues)[1:500][xx] })
Clusters_df <- as.data.frame(Clusters)




#--- Enrichment function per cluster and ontology ---#
run_enrichment <- function(cluster_id, ontology, out_prefix = "GO_enrichment") {
  genes_in_cluster <- cluster_df$Gene[cluster_df$Cluster == cluster_id]
  
  geneList <- factor(as.integer(geneUniverse %in% genes_in_cluster))
  names(geneList) <- geneUniverse
  
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                geneSel = function(p) p == 1,
                annot = annFUN.gene2GO,
                gene2GO = gene2go,
                nodeSize = 10)
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  goResults <- GenTable(GOdata, 
                        classicFisher = resultFisher,
                        topNodes = length(score(resultFisher)))
  
  fname <- paste0(out_prefix, "_Cluster", cluster_id, "_", ontology, ".csv")
  write.csv(goResults, fname, row.names = FALSE)
  
  return(goResults)
}

#--- Loop over clusters and ontologies ---#
results <- list()
for (k in sort(unique(clusters))) {
  for (ont in c("BP", "MF", "CC")) {
    key <- paste0("Cluster", k, "_", ont)
    results[[key]] <- run_enrichment(k, ont)
  }
}