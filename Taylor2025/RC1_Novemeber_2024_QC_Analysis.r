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
> summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   2478   65092  120211  163502  251935  49473

>  summary(colSums(A > 0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1132    9220   10903   10972   13886   15910
####

Getting Reads per UMI Calculation from Summary Files 
reads = read.table('C:/Users/taylo/Desktop/RC1_November2024/Mapped_Data/stringtie_out/read_counts.tab.summary', header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- reads[, !grepl("_unsorted\\.bam$", colnames(reads))]
## colnames(reads) <- sub("^Mapped_Data\\.hisat2_out\\.Justin_", "", colnames(reads))
colnames(reads) <- sub("\\.bam$", "", colnames(reads))

## Getting UMI Counts 
RperU = (reads[1,])/(colSums(A))
write.csv(RperU, "ReadsperUMICalculation_RC1November2024Sequencing.csv")

