#### NTS1 Samples 

LayoutBC = c(rep(NA, 16), 
             rep(c("Tube1", "Tube2"), length.out = 176),  
             rep(c("Tube3", "Tube4"), length.out = 176),  
             rep(NA, 16)) 

sample_labels_HALF1 <- paste(rep(9:96, each = 2), "s", sep = "")  
LayoutBC[17:192] <- paste(LayoutBC[17:192], sample_labels_HALF1, sep = "_")  

sample_labels_HALF2 <- paste(rep(1:88, each = 2), "s", sep = "")   
LayoutBC[193:368] <- paste(LayoutBC[193:368], sample_labels_HALF2, sep = "_") 

convertMatrixtoPlate = function(xx) {
	xx = matrix(xx, ncol = 24)
	rownames(xx) = LETTERS[1:16]
	colnames(xx) = 1:24
	return(xx)
}


LayoutPlate = convertMatrixtoPlate(LayoutBC)


TFInfo = read.csv("NTS1_Geneinfo.csv")


setwd("C:/Users/taylo/Desktop/2025_Nelms/NTS1/")
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

#colnames(A) <- sub("^(\\S+)_.*_(\\d+s)$", "\\1_\\2", colnames(A))

colnames(A) = paste('Tube', sub('[AB]_NTS.+_', '_', sub('.bam.+','',colnames(A))), sep = '')


A <- A[, -16]
colnames(A) 
A[is.na(A)] = 0
A = A[rowSums(A) > 0,]
rownames(A) = annots[rownames(A)]

A <- A[grepl("^Zm", rownames(A)), ]

for (g in unique(rownames(A)[duplicated(rownames(A))])) {
    i = which(rownames(A) == g)
    A[i[1],] = colSums(A[i,])
    A = A[-i[-1],]
}

dim(A)


reord = match(LayoutBC, colnames(A))

library(ComplexHeatmap)
svg('UMI heatmap.svg', width = 8, height =5)
Heatmap(log10(convertMatrixtoPlate(colSums(A,na.rm=T)[reord])), cluster_columns = F, cluster_rows = F)
dev.off()

pseudocount = 100
B2 = A[,colSums(A) >= 5000]
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
B4 = B3[rowSums(B2 >= 10) >= 2,]
logTPM = t(scale(t(B4)))


Heatmap(cor(logTPM))



summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      5    2928   33381   36911   56041  173333
summary(colSums(A > 0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      5    1609    7828    6829   10383   15000



