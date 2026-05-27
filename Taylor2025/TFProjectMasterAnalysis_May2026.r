## TF Project Master Document, R Analysis 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments")

## List of All TFS 
All_GeneName <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/data.csv")

## Maize Gene IDs 

GeneID = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/B73v5_to_B73v4.csv")

## NTS2 metadata 
NTS2_metadata <- read.csv("NTS2_NTS2_1_Metadata.csv")

## NTS3 metadata 

NTS3_metadata <- read.csv("metadata_NTS3_NTS3_1.csv")

## ChipSeq Target Data 
CHIP <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Hu et al/mmc9.csv", header = TRUE, skip = 1, stringsAsFactors = FALSE)
head(CHIP)

## NTS2_UMICounts

A2 <- read.csv("A2_NTS2_NTS2_1.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

NTS2_TFInfo <- read.csv("TFInfo_Mod.csv")

LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

controls = c('A02','P02','B03','O03','H04', 'E06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

## sample_to_well <- setNames(NTS2_metadata$well, NTS2_metadata$sample)
sample_to_well<- setNames(NTS2_metadata$sample, NTS2_metadata$well)
colnames(A2) <- sample_to_well[colnames(A2)]

convertMatrixtoPlate <- function(xx) {
padded <- c(rep(NA, 16), xx, rep(NA, 16))
padded <- padded[1:384]
mat <- matrix(padded, nrow = 16)
rownames(mat) <- LETTERS[1:16]
colnames(mat) <- 1:24
return(mat)
}

A2[is.na(A2)] = 0

A2 <- A2[rowSums(A2) > 0, ]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A2)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

A2 <- A2[, !is.na(colnames(A2)), drop = FALSE]

A2 <- A2[, colSums(!is.na(A2)) > 0]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(NTS2_TFInfo[,9], ', '))),]

sample_to_well <- setNames(NTS2_metadata$well, NTS2_metadata$sample)
colnames(A2) <- sample_to_well[colnames(A2)]

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
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)

library(topGO)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(ComplexHeatmap)

Z = normalize() 
Z_filtered <- Z[rowSums(abs(Z) >= log10(10)) > 0, ]

NTS2_metadata <- NTS2_metadata %>%
  left_join(All_GeneName, by = c("TF" = "gene.ID"))

head(NTS2_metadata)
Samples_HM3 <- c("ZmHB1", "ZmMYB1", "ZmbHLH49", "ZmGATA26", "ZmGLK13", "ZmTHX4")


selected_metadata <- NTS2_metadata %>%
  filter(protein.name %in% Samples_HM3)

selected_wells <- selected_metadata$well
mapped_wells <- LayoutPlate[selected_wells]
wells_to_use <- names(mapped_wells)
replicates_of_wells_to_use <- unname(mapped_wells)
pairs <- t(apply(cbind(wells_to_use, replicates_of_wells_to_use), 1, sort))
unique_pairs <- unique(pairs)
wells_ordered <- as.vector(t(unique_pairs)) 

mCherry2 <- unname(controls)
paired_mCherry2 <- as.vector(rbind(mCherry2[1:12], mCherry2[13:24]))

Z_subset_ordered <- Z[, wells_ordered]
Z_subset_ordered_2 <- Z[, paired_mCherry2]
Z_subset_ordered_3 <- cbind(Z_subset_ordered, Z_subset_ordered_2)

Well_GeneName <- NTS2_metadata %>% 
  filter(well %in% colnames(Z_subset_ordered)) %>%
  select(well, protein.name) %>%
  distinct() %>%           
  deframe()               

col_labels <- Well_GeneName[colnames(Z_subset_ordered_3)]
col_labels[is.na(col_labels)] <- colnames(Z_subset_ordered_3)[is.na(col_labels)]

Z_subset_ordered_4 <- Z_subset_ordered_3[rowSums(abs(Z_subset_ordered_3) >= log10(4)) > 0, ]

min <- min(Z_subset_ordered_4, na.rm = TRUE)
max <- max(Z_subset_ordered_4, na.rm = TRUE)
library(circlize)
hmcols2 = colorRamp2(c(min, 0, max), colors = c('blue', '#EEEEEE', 'red'))

svg("Example_TF_Expression_NTS2_NTS2_1_PlantCenter.svg", width = 15, height = 45)
Heatmap(Z_subset_ordered_4,
        name = "Z-score",
        column_labels = col_labels, 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "TF Examples",
        row_title = "Genes",
        heatmap_legend_param = list(title = "Log10 Fold Change"),
        col = hmcols2)

dev.off()


## NTS3_UMICounts 

NTS3_UMICounts <- read.csv("A2_NTS3_NTS3_1.csv")
NTS3_TFInfo <- read.csv("TFInfo_NTS3_NTS3_1.csv")


LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

controls = c('A02','P02','B03','O03','H04', 'A06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

## sample_to_well <- setNames(NTS3_metadata$well, NTS3_metadata$sample)
sample_to_well<- setNames(NTS3_metadata$sample, NTS3_metadata$well)
colnames(A2) <- sample_to_well[colnames(A2)]

convertMatrixtoPlate <- function(xx) {
padded <- c(rep(NA, 16), xx, rep(NA, 16))
padded <- padded[1:384]
mat <- matrix(padded, nrow = 16)
rownames(mat) <- LETTERS[1:16]
colnames(mat) <- 1:24
return(mat)
}

A2[is.na(A2)] = 0

A2 <- A2[rowSums(A2) > 0, ]
A2 <- A2[grepl("^Zm", rownames(A2)), ]

for (g in unique(rownames(A2)[duplicated(rownames(A2))])) {
    i = which(rownames(A2) == g)
    A2[i[1],] = colSums(A2[i,])
    A2 = A2[-i[-1],]
}

A2 <- A2[, !is.na(colnames(A2)), drop = FALSE]

A2 <- A2[, colSums(!is.na(A2)) > 0]

A2 <- A2[!(rownames(A2) %in% unlist(strsplit(NTS3_TFInfo[,9], ', '))),]

sample_to_well <- setNames(NTS2_metadata$well, NTS2_metadata$sample)
colnames(A2) <- sample_to_well[colnames(A2)]

TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls = c('A02','P02','B03','O03','H04','A06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

Tubes = list(A = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	B = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, c(paste('0', 2:9, sep = ''), 10:12), sep = '') }, simplify=F)),
	C = unlist(sapply(LETTERS[seq(1,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)),
	D = unlist(sapply(LETTERS[seq(2,16,2)], function(xz) { paste(xz, 13:23, sep = '') }, simplify=F)))

normalize = function(crossval = NULL, Z = logTPM, filt = rowSums(A2 >= 10) >= 2, tubespecific = T) {
      if (is.null(crossval)) {
            cnt = controls
      } else {
            cnt = controls[-crossval]
      }

      for (i in 1:length(Tubes)) {  # subtract the mean logTPM value for the controls; this was run separately for each of the four tubes to reduce tube-specific biases
            Z[, colnames(Z) %in% Tubes[[i]]] = sweep(Z[, colnames(Z) %in% Tubes[[i]]], 1, rowMeans(Z[, (colnames(Z) %in% Tubes[[i]]) & (colnames(Z) %in% cnt)]), '-')
      }

      if (is.null(crossval)) {
            return(Z[filt,])
      } else {
            return(Z[filt,controls[!(controls %in% cnt)]])
      }
}

set.seed(1)
CVs2 = sapply(1:choose(24,2), function(i) { normalize(combn(1:24,2)[,i]) }, simplify=F)  # run the cross validation excluding each of the 12 control pairs, one-by-one (leave 1 out cross validation)
corNull2 = unlist(lapply(CVs2, function(xx) { cor(xx[,1], xx[,2]) }))  # correlation between the replicate control samples, using the controls that were not used for normalization during cross validation
sameTubeCV = apply(sub('_.+','',combn(paste(sapply(controls, function(xx) { LETTERS[which(unlist(lapply(Tubes, function(yy) { xx %in% yy })))] }),1:24, sep = '_'),2)), 2, function(xx2) { xx2[1] == xx2[2] })

corNull2[!sameTubeCV]  # << this is what is plotted as the negative controls in that histogram. !sameTubeCV removes control pairs that were both from the same tube, as the normalization is noisy when dropping the number of control samples by 2 (a more sophisticated approach should be less prone to this)


##### Chip Data 
convert_chip_ids <- function(CHIP, GeneID) {
  v5 <- GeneID[,1]              
  v4 <- GeneID[,-1]             
  
  gene_map <- data.frame(
    V4_ID = as.vector(as.matrix(v4)),
    V5_ID = rep(v5, times = ncol(v4)),
    stringsAsFactors = FALSE
  )
  gene_map <- gene_map[ gene_map$V4_ID != "" & !is.na(gene_map$V4_ID), ]
  
  map_list <- split(gene_map$V5_ID, gene_map$V4_ID)
  
  CHIP$V5_regulator <- sapply(CHIP$Regulator, function(x) {
    paste(unique(map_list[[x]]), collapse = ";")
  })
  
  CHIP$V5_target <- sapply(CHIP$Target, function(x) {
    paste(unique(map_list[[x]]), collapse = ";")
  })
  
  return(CHIP)
}

CHIP_converted <- convert_chip_ids(CHIP, GeneID)
head(CHIP_converted)

common_cols <- intersect(names(NTS2_metadata), names(NTS3_metadata))

NTS2_sub <- NTS2_metadata[ , common_cols]
NTS3_sub <- NTS3_metadata[ , common_cols]

TFs <- rbind(NTS2_sub, NTS3_sub)

head(TFs)

IAA10: ZM00001eb156910
IAA27: ZM00001eb336930
IAA29: GRMZM2G141205

Auxin <- c("ZM00001eb156910", "ZM00001eb336930","GRMZM2G141205")
subset_TFs_Auxin <- TFs[ TFs$TF %in% Auxin, ]

common_ids <- intersect(TFs$TF, CHIP_converted$V5_regulator)

subset_chip <- CHIP_converted[ CHIP_converted$V5_regulator %in% common_ids, ]
subset_TFs <- TFs[ TFs$TF %in% common_ids, ]

head(subset_chip)

subset_chip <- CHIP_converted[ CHIP_converted$V5_regulator %in% common_ids, ]
dim(subset_chip)   # should now reflect the expected 14 rows
head(subset_chip)

reg_to_targets_df <- aggregate(V5_target ~ V5_regulator, subset_chip,
                               function(x) unique(x))

dim(reg_to_targets_df)

######################## Forming metadata

library(dplyr)
library(stringr)

TFInfo <- read.csv("TFInfo_Mod.csv")
All_GeneName <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/data.csv")

IDs <- match(TFInfo$Gene.model, All_GeneName$gene.ID)

TFInfo$Gene.model <- All_GeneName$gene.ID[IDs]
TFInfo$GeneName   <- All_GeneName$protein_name[IDs]


TFInfo$Gene.model <- ifelse(
  is.na(TFInfo$Gene.model) & !is.na(TFInfo$V3),
  TFInfo$V3,
  TFInfo$Gene.model
)

LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

plate_rows <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
plate_cols <- sprintf("%02d", 2:23)

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

length(wells) 

barcodes <- rep(paste0(1:88, "s"), times = 4)

meta <- data.frame(
  well = wells,
  BC   = barcodes
)

bc <- paste0(1:88, "s")

AB <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS2_NTS2_1_A_", bc), paste0("NTS2_NTS2_1_B_", bc))
}))

CD <- unlist(lapply(bc, function(bc) {
  c(paste0("NTS2_NTS2_1_C_", bc), paste0("NTS2_NTS2_1_D_", bc))
}))

Wellorder <- c(AB, CD)

A2 <- A2[, Wellorder]


controls = c('A02','P02','B03','O03','H04','E06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

metadata <- data.frame(
  sample = colnames(A2),
  well   = wells
)

metadata$control = metadata$well %in% controls
rownames(metadata) = metadata$sample

TFpos = rep(TFInfo$Gene.model,2)
names(TFpos) = c(TFInfo$Well, LayoutPlate[TFInfo$Well])
metadata$TF = TFpos[metadata$well]

sample_to_well <- setNames(metadata$well, metadata$sample)

colnames(A) <- sample_to_well[colnames(A)]
write.csv(TFInfo, "TFInfo_Mod.csv")

metadata <- metadata %>%
  mutate(batch = str_extract(sample, "(?<=_)[A-Z](?=_[^_]+$)"))
write.csv(metadata, "metadata.csv")