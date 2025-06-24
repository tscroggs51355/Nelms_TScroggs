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

## colnames(A) <- gsub("^Justin_|\\.bam\\.tsv$", "", colnames(A))
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

