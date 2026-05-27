setwd("C:/Users/taylo/Desktop/TF Project Master Doucments") 

metadata = read.csv("Metadata_FullPlate_Jan2026.csv")

A2_NTS2_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS2_NTS2_1.csv", row.names = 1, check.names = FALSE)
A2_NTS3_NTS3_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS3_NTS3_1.csv", row.names = 1, check.names = FALSE)
A2_NTS4_NTS4_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS4_NTS4_1_PreNormalization.csv", row.names = 1, check.names = FALSE)
A2_NTS1_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS1_NTS2_1_PreNormalization.csv", row.names = 1, check.names = FALSE)

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

A2 <- A2_NTS2_NTS2_1[!(rownames(A2_NTS2_NTS2_1) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS2-NTS2"]]


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


Z_NTS2_NTS2_1 = normalize() 

## 
A2 <- A2_NTS3_NTS3_1[!(rownames(A2_NTS3_NTS3_1) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS3-NTS3"]]


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


Z_NTS3_NTS3_1= normalize() 

## 

A2 <- A2_NTS4_NTS4_1[!(rownames(A2_NTS4_NTS4_1) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS4-NTS4"]]


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

Z_NTS4_NTS4_1= normalize() 

##
A2 <- A2_NTS1_NTS2_1[!(rownames(A2_NTS1_NTS2_1) %in% unlist(strsplit(metadata[,5], ', '))),]
TPM = sweep(A2, 2, colSums(A2), '/')*10^6
logTPM = log10(TPM+100)
logTPMf = logTPM[rowSums(A2 >= 10) >= 2,]

controls <- controls_by_plate[["NTS1-NTS2"]]


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

Z_NTS1_NTS2_1= normalize() 

#######
#write.csv(Z_NTS2_NTS2_1, "Z_NTS2_NTS2_1_Normalized_Unfiltered.csv")
#write.csv(Z_NTS3_NTS3_1, "Z_NTS3_NTS3_1_Normalized_Unfiltered.csv")
#write.csv(Z_NTS4_NTS4_1, "Z_NTS4_NTS4_1_Normalized_Unfiltered.csv")
#write.csv(Z_NTS1_NTS2_1, "Z_NTS1_NTS2_1_Normalized_Unfiltered.csv")

Z_NTS2_NTS2_1 = read.csv("Z_NTS2_NTS2_1_Normalized_Unfiltered.csv")
Z_NTS3_NTS3_1 = read.csv("Z_NTS3_NTS3_1_Normalized_Unfiltered.csv")
Z_NTS4_NTS4_1 = read.csv("Z_NTS4_NTS4_1_Normalized_Unfiltered.csv")
Z_NTS1_NTS2_1 = read.csv("Z_NTS1_NTS2_1_Normalized_Unfiltered.csv")


# Start clean
metadata$Quadrant <- NA

# Case 1: hyphen format, e.g. "NTS1-NTS2-1-A_25s"
idx_hyphen <- grepl("^[^-]+-[^-]+-[^-]+-[A-D]_.*$", metadata$SampleName)
metadata$Quadrant[idx_hyphen] <-
  sub("^[^-]+-[^-]+-[^-]+-([A-D])_.*$", "\\1", metadata$SampleName[idx_hyphen])

# Case 2: underscore format, e.g. "NTS2_NTS2_1_A_1s"
idx_underscore <- grepl("^[^_]+_[^_]+_[^_]+_[A-D]_.*$", metadata$SampleName)
metadata$Quadrant[idx_underscore] <-
  sub("^[^_]+_[^_]+_[^_]+_([A-D])_.*$", "\\1", metadata$SampleName[idx_underscore])

# If you already had some clean single-letter quadrants, preserve them:
is_single_letter <- nchar(metadata$Quadrant) == 1 & metadata$Quadrant %in% c("A","B","C","D")
# (If you want to keep any existing ones, you can skip overwriting them in the two blocks above.)

# Sanity check
table(metadata$Quadrant, useNA = "ifany")


############################################################

meta_lookup <- metadata[, c("Assay.Plate", "Array.Well", "SampleID", "Quadrant")]
















Z_NTS2_NTS2_1_annot <- annotate_Z_with_sampleID(Z_NTS2_NTS2_1, "Z_NTS2_NTS2_1", metadata)
Z_NTS3_NTS3_1_annot <- annotate_Z_with_sampleID(Z_NTS3_NTS3_1, "Z_NTS3_NTS3_1", metadata)
Z_NTS4_NTS4_1_annot <- annotate_Z_with_sampleID(Z_NTS4_NTS4_1, "Z_NTS4_NTS4_1", metadata)
Z_NTS1_NTS2_1_annot <- annotate_Z_with_sampleID(Z_NTS1_NTS2_1, "Z_NTS1_NTS2_1", metadata)




datasets <- list(
  NTS2 = Z_NTS2_NTS2_1_annot,
  NTS3 = Z_NTS3_NTS3_1_annot,
  NTS4 = Z_NTS4_NTS4_1_annot
)

reg_counts <- lapply(names(datasets), function(name) {
  df <- datasets[[name]]
  data.frame(
    dataset      = name,
    TF           = colnames(df),
    upregulated  = colSums(df >= log10(2)),
    downregulated = colSums(df <= -log10(2))
  )
})

reg_counts <- do.call(rbind, reg_counts)
reg_counts

reg_counts <- reg_counts[reg_counts$TF != "mCherry", ] ## remove mCherry


library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

reg_counts$total <- reg_counts$upregulated + reg_counts$downregulated
reg_counts$TF <- reorder(reg_counts$TF, -reg_counts$total)

reg_counts_long <- pivot_longer(reg_counts, 
                                 cols = c("upregulated", "downregulated"), 
                                 names_to = "direction", 
                                 values_to = "count")

svg("StackedRegulation.svg", width = 6, height = 6)
ggplot(reg_counts_long, aes(x = TF, y = count, fill = direction)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("upregulated" = "red", "downregulated" = "blue")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "TF", y = "Number of Genes", fill = "Direction")
dev.off()

  #on average, 96 genes upregulated, 67 genes downregulated 

####
metadata_plate = metadata

process_plate <- function(Z, plate_name, metadata_plate) {
  
  plate_label   <- paste0(plate_name, "-", plate_name)
  meta_sub      <- metadata_plate[metadata_plate$Assay.Plate == plate_label, ]
  control_ids   <- meta_sub$SampleID[meta_sub$SampleID == "mCherry"]

  # Split into controls and TFs
  Z_controls <- Z[, colnames(Z) %in% control_ids,  drop = FALSE]
  Z_TF       <- Z[, !colnames(Z) %in% control_ids, drop = FALSE]

  # Strip .1 suffix to get base SampleID for grouping replicates
  sample_ids <- sub("\\.1$", "", colnames(Z_TF))

  # Average the 2 replicates per SampleID using rowMeans
  unique_samples <- unique(sample_ids)
  TF_avg <- sapply(unique_samples, function(s) {
    cols <- which(sample_ids == s)
    rowMeans(Z_TF[, cols, drop = FALSE])
  })
  TF_avg <- as.matrix(TF_avg)
  rownames(TF_avg) <- rownames(Z_TF)
  colnames(TF_avg) <- unique_samples

  # Filter: keep genes upregulated in at least 1 TF
  genes_pass      <- rownames(TF_avg)[rowSums(TF_avg >= log10(2)) > 0]
  TF_avg_filtered <- TF_avg[genes_pass, , drop = FALSE]

  list(
    TF_avg          = TF_avg,
    TF_avg_filtered = TF_avg_filtered,
    controls        = Z_controls,
    genes_pass      = genes_pass,
    meta            = meta_sub
  )
}

processed <- lapply(names(datasets), function(name) {
  process_plate(datasets[[name]], name, metadata_plate)
})
names(processed) <- names(datasets)

# Check dims - should now be 164 TFs per plate
lapply(processed, function(x) dim(x$TF_avg))

# Use filtered genes per plate, then find common ones
common_genes_filtered <- Reduce(intersect, lapply(processed, function(x) x$genes_pass))
cat("Common filtered genes:", length(common_genes_filtered), "\n")

TF_avg_combined <- do.call(cbind, lapply(processed, function(x) {
  x$TF_avg[common_genes_filtered, , drop = FALSE]  # filtered
}))
TF_avg_combined <- as.data.frame(TF_avg_combined)
dim(TF_avg_combined)


# Combine control matrices across all plates using common filtered genes
Z_controls_combined <- do.call(cbind, lapply(processed, function(x) {
  common <- intersect(common_genes_filtered, rownames(x$controls))
  x$controls[common, , drop = FALSE]
}))

# Combine TF and controls side by side
Z_combined_all <- cbind(
  TF_avg_combined[common_genes_filtered, , drop = FALSE],
  Z_controls_combined[common_genes_filtered, , drop = FALSE]
)
Z_combined_all <- as.data.frame(Z_combined_all)


genes_of_interest <- c(
  "Zm00001eb054440",  # TB1
  "Zm00001eb067310",  # WUS1
  "Zm00001eb148390",  # WUS2
  "Zm00001eb144510"   # BBM
)

# Check each plate
lapply(names(processed), function(name) {
  found <- genes_of_interest %in% processed[[name]]$genes_pass
  cat(name, "\n")
  print(data.frame(
    gene    = c("TB1", "WUS1", "WUS2", "BBM"),
    gene_id = genes_of_interest,
    passes  = found
  ))
  cat("\n")
})

** No TB1, WUS, or BBM targets that are above the threshold 

classicalgenes <- read.csv("genes_classical.csv")

# Which classical genes are in TF_avg_combined
matches <- classicalgenes[classicalgenes$v5.Gene.Model.ID %in% rownames(Z_combined_all), 
                          c("v5.Gene.Model.ID", "Gene.Symbol", "Full.Name")]

# Subset to matched genes
Z_classical <- Z_combined_all[matches$v5.Gene.Model.ID, , drop = FALSE]

# Replace gene IDs with gene symbols as row names
gene_symbols <- matches$Gene.Symbol[match(matches$v5.Gene.Model.ID, matches$v5.Gene.Model.ID)]
rownames(Z_classical) <- make.unique(gene_symbols)

# Scale by log10(2)
Z_classical_scaled <- Z_classical / log10(2)

# Define column groups - TF columns vs Control columns
column_group <- c(
  rep("TF", ncol(TF_avg_combined)),
  rep("Control", ncol(Z_controls_combined))
)

# Plot heatmap
hmcols <- colorRamp2(c(-4, 0, 4), colors = c('blue', '#EEEEEE', 'red'))

svg("ClassicalExpression.svg", width = 5, height = 10)
Heatmap(
  as.matrix(Z_classical_scaled),
  name = "Expression",
  col = hmcols,
  column_split = column_group,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_gap = unit(5, "mm"),
  row_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  )
)
dev.off()