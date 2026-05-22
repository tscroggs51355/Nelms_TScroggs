# ============================================================
# August 2024 Maize RNA-seq Analysis: R, C1, R+C1 vs Control
# ============================================================

setwd("C:/Users/taylo/Desktop/AugustSequencingNelms2024/")

# ── 1. LOAD GENE ANNOTATIONS ────────────────────────────────
# Parse the merged StringTie GTF to build a gene_id → gene_name lookup.
# Each line's 9th field contains semicolon-delimited attributes; we split those
# and pull out the last attribute (gene name when 3 fields are present,
# otherwise the gene_id itself). Duplicates are dropped so names stay unique.

annots <- strsplit(read.table('Mapped_Data/stringtie_out/stringtie_merged.gtf',
                               sep = '\t')[, 9], '; ')
names(annots) <- unlist(lapply(annots, function(xx) xx[1]))
names(annots) <- sub('gene_id ', '', names(annots))
annots <- annots[!duplicated(names(annots))]
annots <- sub(';', '', sub(' ', '',
               unlist(lapply(annots, function(xx) {
                 sub('.+ ', '', if (length(xx) == 3) xx[3] else xx[1])
               }))))


# ── 2. READ UMI COUNT FILES ──────────────────────────────────
# Each .tsv in Mapped_Data/UMIcounts/ is one sample (one BAM file).
# Files are loaded into a list, then merged into a single count matrix A
# where rows = genes and columns = samples. Missing values become 0.

files <- dir('Mapped_Data/UMIcounts')

A0 <- list()
for (f in files) {
  A0[[f]] <- read.table(paste0('Mapped_Data/UMIcounts/', f),
                         sep = '\t', header = TRUE, row.names = 1)
}

# Union of all gene IDs across samples
gn <- unique(unlist(lapply(A0, rownames)))

A <- matrix(NA, nrow = length(gn), ncol = length(files))
rownames(A) <- gn
colnames(A) <- files

for (f in files) {
  A[, f] <- A0[[f]][match(gn, rownames(A0[[f]])), 1]
}


# ── 3. CLEAN COLUMN NAMES AND MATRIX ────────────────────────
# Strip the filename down to just the sample number + 's' (e.g. "23s").
# Replace NA with 0, drop all-zero gene rows, and swap gene IDs for names.

colnames(A) <- sub(".*_(\\d+s)\\.bam\\.tsv$", "\\1", colnames(A))
A[is.na(A)] <- 0
A <- A[rowSums(A) > 0, ]
rownames(A) <- annots[rownames(A)]


# ── 4. COLLAPSE DUPLICATE GENE NAMES ────────────────────────
# StringTie can assign the same name to multiple loci. Where that happens,
# sum the counts across duplicates into a single row then remove the extras.

for (g in unique(rownames(A)[duplicated(rownames(A))])) {
  i <- which(rownames(A) == g)
  A[i[1], ] <- colSums(A[i, ])
  A <- A[-i[-1], ]
}


# ── 5. FILTER TO MAIZE GENES AND REMOVE BAD SAMPLE ──────────
# Keep only annotated Zea mays genes (prefix "Zm").
# Sample 19s is flagged as a quality outlier and dropped before reordering.

A <- A[grepl("^Zm", rownames(A)), ]   # 26,257 genes × 96 samples
A <- A[, -which(colnames(A) == "19s")] # drop outlier sample 19s
A <- A[, order(as.numeric(sub('s', '', colnames(A))))]  # sort columns numerically
A <- A[, 1:31]   # retain only samples 1–31 (drop 32–96; unused groups/lanes)


# ── 6. DEFINE SAMPLE GROUPS AND BUILD DESEQ2 INPUT ──────────
# Assign each column to its experimental condition (treatment group).
# Reorder the count matrix to match the sample order used in colData.

samples <- c(
  "5s",  "6s",  "7s",  "8s",  "25s", "26s", "27s", "28s",  # R        (n=8)
  "17s", "18s", "20s", "21s", "22s", "23s", "24s",           # C1       (n=7)
  "1s",  "2s",  "3s",  "4s",  "9s",  "10s", "11s", "12s",   # R+C1     (n=8)
  "13s", "14s", "15s", "16s", "29s", "30s", "31s"            # Control  (n=7 after 19s removal)
)

conditions <- c(rep("R",       8),
                rep("C1",      7),
                rep("R+C1",    8),
                rep("Control", 7))

colData <- data.frame(
  condition = factor(conditions, levels = c("R", "C1", "R+C1", "Control"))
)
rownames(colData) <- samples

counts <- A[, samples]   # reorder columns to match colData
counts[is.na(counts)] <- 0


# ── 7. NORMALIZED EXPRESSION MATRIX FOR CLUSTERING ──────────
# B2: samples with ≥5,000 total counts (library-size filter).
# B3: log10(TPM + 100) — pseudocount of 100 reduces noise from low-count genes.
# B4: restrict to genes with ≥10 counts in at least 2 samples.
# logTPM: z-score each gene across samples (mean=0, sd=1) for clustering.

pseudocount <- 100
B2 <- A[, colSums(A) >= 5000]
B3 <- log(sweep(B2, 2, colSums(B2), '/') * 1e6 + pseudocount, 10)
B4 <- B3[rowSums(B2 >= 10) >= 2, ]   # genes expressed above threshold in ≥2 samples
logTPM <- t(scale(t(B4)))             # z-score across samples per gene


# ── 8. DIFFERENTIAL EXPRESSION WITH DESEQ2 ──────────────────
# Fit a negative binomial model with condition as the sole factor.
# Three pairwise contrasts are computed, each comparing a treatment to Control.

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = colData,
                               design    = ~ condition)
dds <- DESeq(dds)

R_control     <- results(dds, contrast = c("condition", "R",     "Control"))
C1_control    <- results(dds, contrast = c("condition", "C1",    "Control"))
R_C1_control  <- results(dds, contrast = c("condition", "R+C1",  "Control"))


# ── 9. P-VALUE POOLING (FISHER'S METHOD) ────────────────────
# Combine per-gene p-values across the three contrasts using Fisher's
# chi-squared statistic: X² = −2 · Σ ln(p).
# Degrees of freedom = 2 × number of tests (6 here).
# NOTE: confirm whether adding dchisq() is appropriate before publication
#       (Brad to review this adjustment term).

combined_pvalues <- data.frame(
  R    = R_control$pvalue,
  C1   = C1_control$pvalue,
  R_C1 = R_C1_control$pvalue,
  row.names = rownames(R_control)
)

X2 <- -2 * rowSums(log(combined_pvalues))
combined_pvalues$p_pool <- pchisq(X2, df = 2 * ncol(combined_pvalues), lower.tail = FALSE) +
                            dchisq(X2, df = 2 * 3)   # ← Brad to verify this term


# ── 10. MULTIPLE TESTING CORRECTION ─────────────────────────
# Holm correction is used throughout (more conservative than BH/FDR).
# maxLog2: for each gene, pick the log2FC with the largest absolute value
#          across the three contrasts — used as a summary effect size.

combined_pvalues$padj <- p.adjust(combined_pvalues$p_pool, method = "holm")

combined_pvalues$maxLog2 <- apply(
  cbind(R_control[, 2], C1_control[, 2], R_C1_control[, 2]),
  1,
  function(xx) xx[rank(-abs(xx), ties.method = 'first') == 1]
)

# Significant = Holm-adjusted p ≤ 0.05 AND |maxLog2FC| ≥ 1 (≥ 2-fold change)
combined_pvalues$Significant <- (combined_pvalues$padj <= 0.05) &
                                 (abs(combined_pvalues$maxLog2) >= 1)
combined_pvalues$Significant[is.na(combined_pvalues$Significant)] <- FALSE

# Sort: significant genes first, then by ascending adjusted p-value
combined_pvalues <- combined_pvalues[order(-combined_pvalues$Significant,
                                            combined_pvalues$padj), ]
sigGenes <- rownames(combined_pvalues)[combined_pvalues$Significant]


# ── 11. HIERARCHICAL CLUSTERING ─────────────────────────────
# Cluster genes using z-scored logTPM values (NOT raw counts).
# Ward's D2 linkage on Euclidean distances; 4 clusters chosen after
# visual inspection (5 was also considered but 4 retained for GO analysis).
# Only DEGs (sigGenes) are written to the output CSV.

d  <- dist(logTPM, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

plot(hc, labels = rownames(logTPM),
     main = "Hierarchical Clustering Dendrogram", xlab = "Gene Names")

num_clusters <- 4
clusters     <- cutree(hc, k = num_clusters)
# table(clusters): 1=2126, 2=2085, 3=4287, 4=1766

# Subset to significant DEGs and export for GO term analysis
filtered_clusters <- clusters[names(clusters) %in% sigGenes]
clusters_df <- data.frame(gene = names(filtered_clusters),
                           cluster = as.numeric(filtered_clusters))
clusters_df <- clusters_df[order(clusters_df$cluster), ]
write.csv(clusters_df, "clusters_df.csv", row.names = FALSE)