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
## Run 6/5/2025 1:13 to here - 























pseudocount = 100
B2 = A[,colSums(A) >= 5000]
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,]
logTPM = t(scale(t(B4)))
131/136 passed QC 

SubsetColnames <- setdiff(colnames(A), colnames(B4))
print(SubsetColnames)
[1] "T384-1_34s" "T384-1_40s" "T384-1_6s"  "T384-2_15s" "T384-2_25s"


#########################################################
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    961   23208   39590   46931   59655  383608

 summary(colSums(A > 0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    583    3918    4469    4689    5513   15875

#########################################################


library(ComplexHeatmap)

CorrelationMatrixData <- B4[, !colnames(B4) %in% c("T96-1_1s", "T384-1_2s", "T96-1_2s")]
cor_matrix = cor(CorrelationMatrixData, method = "pearson")
Heatmap(cor_matrix)
svg("CorrelationHeatmap_TS_March2025.svg", width = 24, height = 24)
Heatmap(cor_matrix)
dev.off()

#Removed the following after plotting initially, because they were completeley driving the correlation 
#T96-1_1S, T384-1_2s, T96-1_2s 



####

Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/June2025_Sequencing/Mapped_Data/stringtie_out/read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- reads[, !grepl("_unsorted\\.bam$", colnames(reads))]
colnames(reads) <- sub("^Mapped_Data\\.hisat2_out\\.Justin_", "", colnames(reads))
colnames(reads) <- sub("\\.bam$", "", colnames(reads))




## Getting UMI Counts 




RperU = (reads[1,])/(colSums(A))


write.csv(RperU, "ReadsperUMICalculation_June2025Sequencing.csv")

