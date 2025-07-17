## Location of Files 
 /scratch/tms51355/Taylor2025/July2025Sequencing_NTS2


 setwd("C:/Users/taylo/Desktop/2025_Nelms/July2025/NTS2_Seq")


NTS2_Array_Clean <- read.csv("NTS2_Array_clean.csv")
NTS2_Array_Clean$Well <- paste0(NTS2_Array_Clean$Row, sprintf("%02d", NTS2_Array_Clean$Column))
library(dplyr)


NTS2_Array <- NTS2_Array_Clean %>%
  left_join(Batch_Synbio_Combined %>% select(SampleID, TFomeStockID), by = "SampleID")

TFome_Data <- read.csv("Maize_TFome_Bulk_data _ corrected BDN.csv", header = TRUE)
TFome_Data <- TFome_Data[,1:10]

Batch1_Synbio <- read.csv("Batch1_Synbio_NTS2.csv")

Batch2_Synbio <- read.csv("Batch2_Synbio_NTS2.csv")

Batch3_Synbio <- read.csv("Batch3_Synbio_NTS2.csv")
Batch_Synbio_Combined <- bind_rows(Batch1_Synbio, Batch2_Synbio, Batch3_Synbio)
NTS2_Array <- NTS2_Array %>%
  left_join(TFome_Data %>% select(Stock.number, Gene.model),
            by = c("TFomeStockID" = "Stock.number"))
write.csv(NTS2_Array, "NTS2_Geneinfo.csv")

MetaData: 
TFInfo = read.csv("NTS2_Geneinfo.csv") 


LayoutBC = c(rep(c("Tube1", "Tube2"), length.out = 192),  
             rep(c("Tube3", "Tube4"), length.out = 192)) 

sample_labels_HALF1 <- paste(rep(1:96, each = 2), "s", sep = "")  
LayoutBC[1:192] <- paste(LayoutBC[1:192], sample_labels_HALF1, sep = "_")  

sample_labels_HALF2 <- paste(rep(1:96, each = 2), "s", sep = "")   
LayoutBC[193:384] <- paste(LayoutBC[193:384], sample_labels_HALF2, sep = "_") 
head(LayoutBC)


LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]



convertMatrixtoPlate = function(xx) {
	xx = matrix(xx, ncol = 24)
	rownames(xx) = LETTERS[1:16]
	colnames(xx) = 1:24
	return(xx)
}


rows = LETTERS[1:16]
columns = sprintf("%02d", 1:24)
rows_columns = expand.grid(rows,columns)
rows_columns = paste(rows_columns$Var1, rows_columns$Var2, sep= "")



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

A <- A[, -16]
colnames(A) = paste('Tube', sub('[AB]_NTS.+_', '_', sub('.bam.+','',colnames(A))), sep = '')

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