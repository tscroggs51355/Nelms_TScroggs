## August 20th, Remapped Data 

setwd("C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August")


Augreads = read.table('C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August/Mapped_Data/Reads/all_read_counts.tab.summary',  header=T, sep = '\t', stringsAsFactors=F, row.names=1)
Julyreads = read.table('C:/Users/taylo/Desktop/2025_Nelms/August 2025/NTS2_NTS2_1_Sequencing_July_August/Mapped_Data/JulyReads/read_counts.tab.summary',  header=T, sep = '\t', stringsAsFactors=F, row.names=1)

reads <- cbind(Augreads,Julyreads)
write.csv(reads,"reads_updatedmapping.csv")

colnames(reads) <- substr(colnames(reads), 78, nchar(colnames(reads)) - 4)

colnames(reads) <- gsub("_S[0-9]+_L[0-9]+", "", colnames(reads))


colnames(reads)[1:192] <- sapply(colnames(reads)[1:192], function(x) substr(x, 9, nchar(x)))

reads_XX= reads

TFInfo = read.csv("C:/Users/taylo/Desktop/2025_Nelms/July2025/NTS2_Seq/TFInfo_Mod.csv")

metadata = read.csv("C:/Users/taylo/Desktop/2025_Nelms/July2025/NTS2_Seq/metadata.csv")

colnames(reads_XX) = metadata$well[match(colnames(reads_XX), metadata$sample)]
#TFInfo <- merge(TFInfo, metadata[, c("Well", "sample")], by = "Well", all.x = TRUE)
metadata$GeneName <- TFInfo$GeneName[match(metadata$TF, TFInfo$Gene.model)]

## remove 89 - 96 barcodes 
reads_XX <- reads_XX[, !is.na(colnames(reads_XX))]



files = dir('Mapped_Data/UMIcounts')
files <- files[!grepl("_8[9]s\\.tsv$|_9[0-6]s\\.tsv$", files)]

A = list()
for (f in files) {
	A[[f]] = read.table(paste('Mapped_Data/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gns = unique(unlist(lapply(A, rownames)))
A2 = matrix(NA, nrow = length(gns), ncol = length(A))
rownames(A2) = gns
colnames(A2) = files

colnames(A2)[89:264] <- gsub("_S[0-9]+_L[0-9]+", "", colnames(A2)[89:264])


colnames(A2) <- sub("\\.tsv$", "", colnames(A2))

for (i in 1:length(A)) {
    A2[rownames(A[[i]]),i] = A[[i]][,1]
    cat(i,'\n')
}
A2[is.na(A2)] = 0
A2 = A2[rowSums(A2) > 0,]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}



### Getting Into Normalization #### 

convertMatrixtoPlate <- function(xx) {
  padded <- c(rep(NA, 16), xx, rep(NA, 16))
  padded <- padded[1:384]  
  
  mat <- matrix(padded, nrow = 16)
  rownames(mat) <- LETTERS[1:16]
  colnames(mat) <- 1:24

  # Explicitly set I09 to NA
  mat["I", "9"] <- NA

  return(mat)
}


LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

sample_to_well <- setNames(metadata$well, metadata$sample)
colnames(A2) <- sample_to_well[colnames(A2)]
A2 <- A2[,match(LayoutPlate,colnames(A2))]

## Getting UMI Counts 

A2 <- A2[, !is.na(colnames(A2)), drop = FALSE]
valid_cols <- intersect(colnames(A2), colnames(reads_XX))
A2 <- A2[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, valid_cols, drop = FALSE]
reads_XX <- reads_XX[, match(colnames(A2), colnames(reads_XX))]
RperU = (reads_XX[1,])/(colSums(A2))
RperU <- as.data.frame(RperU)
RperU_vec <- as.numeric(RperU["Assigned", ])
names(RperU_vec) <- colnames(RperU)
svg("Brad/ReadPerUMI_NTS2_NTS2_1_Barplot.svg", width = 10, height =12)
barplot(RperU_vec,
        las = 2,              
        col = "hotpink",
        border = NA,
        main = "Reads per UMI",
        ylab = "Reads per UMI",
        cex.names = 0.7, 
        space= 5.0)
dev.off()
write.csv(RperU, "Brad/ReadsPerUMI_V5_NTS2_NTS2_1.csv")



library(ComplexHeatmap)

convertMatrixtoPlate <- function(xx) {
  mat <- matrix(xx, nrow = 16, byrow = FALSE)
  rownames(mat) <- LETTERS[1:16]
  colnames(mat) <- 1:24

   mat["I", "9"] <- NA

  return(mat)
}
A2 <- A2[, colSums(!is.na(A2)) > 0]

summary(colSums(A2))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  48657  222834  346264  406377  546457 2258166 
summary(colSums(A2>0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  10098   14898   16927   16827   18584   23108

convertMatrixtoPlate <- function(xx) {
padded <- c(rep(NA, 16), xx, rep(NA, 16))
padded <- padded[1:384]
mat <- matrix(padded, nrow = 16)
rownames(mat) <- LETTERS[1:16]
colnames(mat) <- 1:24
return(mat)
}
svg("UMICounts_Heatmap.svg", width = 25, height = 25)
Heatmap((log10(convertMatrixtoPlate(colSums(A2,na.rm=T)))), cluster_columns = F, cluster_rows = F)
dev.off()



##write.csv(A2, "Brad/A2_NTS2_NTS2_1.csv")
A2 <- A2[!(rownames(A2) %in% unlist(strsplit(TFInfo[,9], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls = c('A02','P02','B03','O03','H04','E06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
	if (is.null(crossval)) {
			cnt = controls
		} else if (length(crossval) == 1) {
			cnt = controls[-(c(0,12) + crossval)]  # if crossval is not null, remove a single control well and its corresponding replicate across the plate
		} else {
			cnt = controls[-crossval]
		}
	
	if (tubespecific) { 
		for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
			Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
		}
	} else { Z = sweep(Z, 1, rowMeans(Z), '-') }
	
	if (is.null(crossval)) {
			return(Z[filt,])
		} else {
			return(Z[filt,controls[!(controls %in% cnt)]])
		}
}

Z = normalize(tubespecific = T)
cors = cor(Z[rank(-apply(Z,1,sd)) <= 2000,])
diag(cors) = NA 

library(ComplexHeatmap)
library(seriation)
minmax = function(x) {
	sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}
minmax2 = function(x) {
	x = sweep(x, 1, apply(x, 1, min), '-')
	sweep(x, 1, apply(x, 1, max), '/')
}


SScore = function(A, B, alpha = .2, len = NULL) {
	# xx = match(A, B)  ## For named vectores
	xx = rank(-A)[order(-B)]
	#xx = xx[!is.na(A) & !is.na(B)]
	#if (is.null(leng)) { len = max(c(round(9.21034/alpha + 1), length(xx))) }
	xx = sapply(1:len, function(i) { sum(xx[1:i] <= i) })
	return(sum(exp(-alpha*1:length(xx))*xx))
}

SSmat = function(x, alpha = .2) {
            ncx <- ncol(x)
            len = max(c(round(9.21034/alpha + 1), nrow(x)))
            r <- matrix(0, nrow = ncx, ncol = ncx)
            for (i in seq_len(ncx)) {
                for (j in i:ncx) {
		  x2 = x[,i]
		  y2 = x[,j]
                  ok <- complete.cases(x2, y2)
                  r[i, j] <- if (any(ok)) {
                    SScore(x2[ok], y2[ok], alpha = alpha, len = len)
                  } else NA
                  r[j, i] = r[i, j]
                }
            }
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
	return(r)
}


hmcols = colorRampPalette(c('#0571b0','#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020','#ca0020'))(100)
o1 = seriate(dist(cors^3), method = "OLO")

## Plotting correlatin of all samples
svg("Brad/Correlation_AllSamples_NTS2_NTS2_1.svg", width = 30, height = 30)
hm=Heatmap(cors, col = hmcols, use_raster = T, raster_device = 'png', show_heatmap_legend = T, show_row_names = T, show_column_names = T, cluster_columns = as.dendrogram(o1[[1]]), cluster_rows = as.dendrogram(o1[[1]]))
dev.off()

## Plotting correlation of just the controls 
svg("Brad/Correlation_Controls_NTS2_NTS2_1.svg", width = 30, height = 30)
control_label <- unique(as.character(controls))
hm= Heatmap(cors[control_label, control_label], col = hmcols, use_raster = T, raster_device = 'png', show_heatmap_legend = T, show_row_names = T, show_column_names = T, cluster_columns = as.dendrogram(o1[[1]]), cluster_rows = as.dendrogram(o1[[1]]))
dev.off()
###############################################################
## Plotting Correlation beween the batches 
##Batch A Correlation Heatmap 

batch_A <- metadata$well[metadata$batch == "A"]
batch_A_samples <- intersect(batch_A, rownames(cors))
gene_names_A <- setNames(TFInfo$GeneName, TFInfo$Well)
labels <- gene_names_A[batch_A_samples]
labels[is.na(labels)] <- batch_A_samples[is.na(labels)]
labels <- make.unique(labels)
cors_batch_A <- cors[batch_A_samples, batch_A_samples]
rownames(cors_batch_A) <- labels
colnames(cors_batch_A) <- labels
hm1 =Heatmap(cors_batch_A, col = hmcols, use_raster = TRUE)


##Batch B Correlation Heatmap 

batch_B <- metadata$well[metadata$batch == "B"]
batch_B_samples <- intersect(batch_B, rownames(cors))
gene_names_B <- setNames(TFInfo$GeneName, TFInfo$Well)
labels <- gene_names_B[batch_B_samples]
labels[is.na(labels)] <- batch_B_samples[is.na(labels)]
labels <- make.unique(labels)
cors_batch_B <- cors[batch_B_samples, batch_B_samples]
rownames(cors_batch_B) <- labels
colnames(cors_batch_B) <- labels
hm2 =Heatmap(cors_batch_B, col = hmcols, use_raster = TRUE)


##Batch C Correlation Heatmap 

batch_C <- metadata$well[metadata$batch == "C"]
batch_C_samples <- intersect(batch_C, rownames(cors))
gene_names_C <- setNames(metadata$GeneName, metadata$well)
labels <- gene_names_C[batch_C_samples]
labels[is.na(labels)] <- batch_C_samples[is.na(labels)]
labels <- make.unique(labels)
cors_batch_C <- cors[batch_C_samples, batch_C_samples]
rownames(cors_batch_C) <- labels
colnames(cors_batch_C) <- labels
hm3 =Heatmap(cors_batch_C, col = hmcols, use_raster = TRUE)

##Batch D Correlation Heatmap 

batch_D <- metadata$well[metadata$batch == "D"]
batch_D_samples <- intersect(batch_D, rownames(cors))
gene_names_D <- setNames(metadata$GeneName, metadata$well)
labels <- gene_names_C[batch_D_samples]
labels[is.na(labels)] <- batch_D_samples[is.na(labels)]
labels <- make.unique(labels)
cors_batch_D <- cors[batch_D_samples, batch_D_samples]
rownames(cors_batch_D) <- labels
colnames(cors_batch_D) <- labels
hm4 = Heatmap(cors_batch_D, col = hmcols, use_raster = TRUE)

library(gridExtra)
library(grid)

g1 <- grid.grabExpr(draw(hm1))
g2 <- grid.grabExpr(draw(hm2))
g3 <- grid.grabExpr(draw(hm3))
g4 <- grid.grabExpr(draw(hm4))

svg("Brad/correlation_batches.svg", width = 25, height = 25)
grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()

##################################################
# Colors

batch_levels <- unique(metadata$batch)
batch_palette <- setNames(
  RColorBrewer::brewer.pal(length(batch_levels), "Set3")[seq_along(batch_levels)],
  batch_levels
)

batch <- setNames(metadata$batch, metadata$well)

top_annot <- ComplexHeatmap::HeatmapAnnotation(
  Batch = batch[rownames(cors)],  # or rownames(cors_batch) etc.
  col = list(Batch = batch_palette),
  annotation_name_side = "left"
)

hm <- Heatmap(
  cors,
  col = hmcols,
  use_raster = TRUE, raster_device = 'png',
  show_heatmap_legend = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,  
  cluster_columns = as.dendrogram(o1[[1]]),
  cluster_rows = as.dendrogram(o1[[1]]),
  top_annotation = top_annot
)

hm_drawn <- draw(hm)

row_labels_ordered <- rownames(cors)[row_order(hm_drawn)]
col_labels_ordered <- colnames(cors)[column_order(hm_drawn)]

hm_final <- Heatmap(
  cors,
  col = hmcols,
  use_raster = TRUE, raster_device = 'png',
  show_heatmap_legend = TRUE,
  show_row_names = TRUE,
  row_labels = row_labels_ordered,
  column_labels = col_labels_ordered,
  show_column_names = TRUE,  
  cluster_columns = as.dendrogram(o1[[1]]),
  cluster_rows = as.dendrogram(o1[[1]]),
  top_annotation = top_annot
)

draw(hm_final)

svg('Brad/correlationheatmap_withannotation.svg', width = 25, height = 25)
draw(hm_final)
dev.off()

















convertMatrixtoPlate <- function(xx) {
  padded <- c(rep(NA, 16), xx, rep(NA, 16))
  padded <- padded[1:384]  
  
  mat <- matrix(padded, nrow = 16)
  rownames(mat) <- LETTERS[1:16]
  colnames(mat) <- 1:24
  return(mat)
}

Heatmap((convertMatrixtoPlate(log10((colSums(A2))))), cluster_columns = F, cluster_rows = F)

svg("ReadsperUMI.svg", width = 40,  height = 10) 
barplot(RperU_df$RperU_matrix, 
        names.arg = RperU_df$Sample, 
        col = "pink", 
        main = "Barplot of RperU_matrix values for each Sample", 
        xlab = "Sample", 
        ylab = "ReadsPerUMI", 
        las = 2)
        
dev.off()

                   