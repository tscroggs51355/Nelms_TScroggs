setwd("C:/Users/taylo/Desktop/TF Project Master Doucments") 

metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Metadata_FullPlate_Jan2026.csv")

### 1. Extract quadrant from SampleName
metadata$Quadrant <- sub(".*-([A-D])_.*", "\\1", metadata$SampleName)

### 2. Build metadata lookup
meta_lookup <- metadata[, c("Assay.Plate", "Array.Well", "SampleID", "SampleName", "Quadrant")]

A2_NTS2_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS2_NTS2_1.csv", row.names = 1, check.names = FALSE)
A2_NTS3_NTS3_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS3_NTS3_1.csv", row.names = 1, check.names = FALSE)
A2_NTS4_NTS4_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS4_NTS4_1_PreNormalization.csv", row.names = 1, check.names = FALSE)
A2_NTS1_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS1_NTS2_1_PreNormalization.csv", row.names = 1, check.names = FALSE)

A2_list <- list(
  NTS2_NTS2_1 = A2_NTS2_NTS2_1,
  NTS3_NTS3_1 = A2_NTS3_NTS3_1,
  NTS4_NTS4_1 = A2_NTS4_NTS4_1
)

Z_NTS2_NTS2_1 = read.csv("Z_NTS2_NTS2_1_Normalized_Unfiltered.csv", row.names = 1, check.names = FALSE)
Z_NTS3_NTS3_1 = read.csv("Z_NTS3_NTS3_1_Normalized_Unfiltered.csv", row.names = 1, check.names = FALSE)
Z_NTS4_NTS4_1 = read.csv("Z_NTS4_NTS4_1_Normalized_Unfiltered.csv", row.names = 1, check.names = FALSE)

rename_to_samplename <- function(df, plate_name, metadata) {

  # Metadata for this plate
  meta_sub <- metadata[metadata$Assay.Plate == plate_name, ]

  # Create lookup: Array.Well -> SampleName
  lookup <- setNames(meta_sub$SampleID, meta_sub$Array.Well)

  # Current column names
  old_names <- colnames(df)

  # Replace wells with SampleName when available
  new_names <- ifelse(
    old_names %in% names(lookup) & !is.na(lookup[old_names]),
    lookup[old_names],
    old_names
  )

  colnames(df) <- new_names
  df
}

z_dfs <- list(
  "NTS2-NTS2" = Z_NTS2_NTS2_1,
  "NTS3-NTS3" = Z_NTS3_NTS3_1,
  "NTS4-NTS4" = Z_NTS4_NTS4_1
)

z_dfs <- lapply(names(z_dfs), function(plate) {
  rename_to_samplename(z_dfs[[plate]], plate, metadata)
})

names(z_dfs) <- c("NTS2-NTS2", "NTS3-NTS3", "NTS4-NTS4")

# Put them back into your workspace
Z_NTS2_NTS2_1 <- z_dfs[["NTS2-NTS2"]]
Z_NTS3_NTS3_1 <- z_dfs[["NTS3-NTS3"]]
Z_NTS4_NTS4_1 <- z_dfs[["NTS4-NTS4"]]

mc_dfs <- lapply(z_dfs, function(df) df["mCherry"])

names(mc_dfs)
# "NTS2-NTS2" "NTS3-NTS3" "NTS4-NTS4"

Z_NTS2_NTS2_1_MC <- mc_dfs[["NTS2-NTS2"]]
Z_NTS3_NTS3_1_MC <- mc_dfs[["NTS3-NTS3"]]
Z_NTS4_NTS4_1_MC <- mc_dfs[["NTS4-NTS4"]]

Z_NTS2_NTS2_1 <- Z_NTS2_NTS2_1[, colnames(Z_NTS2_NTS2_1) != "mCherry"]
Z_NTS3_NTS3_1 <- Z_NTS3_NTS3_1[, colnames(Z_NTS3_NTS3_1) != "mCherry"]
Z_NTS4_NTS4_1 <- Z_NTS4_NTS4_1[, colnames(Z_NTS4_NTS4_1) != "mCherry"]

reorder_pairs <- function(df) {

  cols <- colnames(df)

  # base name (remove .1 if present)
  base <- sub("\\.1$", "", cols)

  # split into paired structure
  ordered_cols <- unlist(lapply(unique(base), function(b) {
    c(
      cols[cols == b],
      cols[cols == paste0(b, ".1")]
    )
  }))

  # keep only existing (in case some pairs missing)
  ordered_cols <- ordered_cols[ordered_cols %in% cols]

  df[, ordered_cols, drop = FALSE]
}

Z_NTS2_NTS2_1 <- reorder_pairs(Z_NTS2_NTS2_1)
Z_NTS3_NTS3_1 <- reorder_pairs(Z_NTS3_NTS3_1)
Z_NTS4_NTS4_1 <- reorder_pairs(Z_NTS4_NTS4_1)

cutoff <- log10(2)

dataset <- list(
  Z_NTS2_NTS2_1, 
  Z_NTS3_NTS3_1, 
  Z_NTS4_NTS4_1
)

upregulated <- lapply(dataset, function(x) {
  apply(x, 2, function(col) {
    rownames(x)[col >= cutoff]
  })
})

n_targets_list <- lapply(seq_along(upregulated), function(i) {
  data.frame(
    Sample = names(upregulated[[i]]),
    N = sapply(upregulated[[i]], length)
  )
})

library(ggplot2)

p1 <- ggplot(n_targets_list[[1]], aes(x = Sample, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("NTS2-NTS2-1")

p2 <- ggplot(n_targets_list[[2]], aes(x = Sample, y = N)) +
  geom_bar(stat = "identity", fill = "orange") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("NTS3-NTS3-1")

p3 <- ggplot(n_targets_list[[3]], aes(x = Sample, y = N)) +
  geom_bar(stat = "identity", fill = "pink") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("NTS4-NTS4-1")


ggsave("NTS2_NTS2_1_upregulated.pdf", p1, width = 44, height = 12)
ggsave("NTS3_NTS3_1_upregulated.pdf", p2, width = 44, height = 12)
ggsave("NTS4_NTS4_1_upregulated.pdf", p3, width = 44, height = 12)


library(grid)

plot_correlation_heatmap <- function(Z, name, out_file) {

  Z_sub <- Z[rank(-apply(Z, 1, sd, na.rm = TRUE)) <= 2000, ]

  cors <- cor(Z_sub, use = "pairwise.complete.obs")
  diag(cors) <- NA

  o1 <- seriate(dist(cors^3), method = "OLO")

  hmcols <- colorRampPalette(
    c('#0571b0','#0571b0','#92c5de','#f7f7f7',
      '#f4a582','#ca0020','#ca0020')
  )(100)

  pdf(out_file, width = 52, height = 52)

  ht <- Heatmap(
    cors,
    col = hmcols,
    use_raster = TRUE,
    raster_device = "png",
    show_heatmap_legend = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_columns = as.dendrogram(o1[[1]]),
    cluster_rows = as.dendrogram(o1[[1]])
  )

  draw(ht)

  dev.off()
}

plot_correlation_heatmap(Z_NTS2_NTS2_1, "NTS2", "Correlation_NTS2.pdf")
plot_correlation_heatmap(Z_NTS3_NTS3_1, "NTS3", "Correlation_NTS3.pdf")
plot_correlation_heatmap(Z_NTS4_NTS4_1, "NTS4", "Correlation_NTS4.pdf")

## 
get_cor_matrix <- function(Z) {
  Z_sub <- Z[rank(-apply(Z, 1, sd, na.rm = TRUE)) <= 2000, ]
  cor(Z_sub, use = "pairwise.complete.obs")
}

get_replicate_cor <- function(cors) {

  base <- sub("\\.1$", "", colnames(cors))

  idx <- which(outer(colnames(cors), colnames(cors),
                     Vectorize(function(a, b)
                       sub("\\.1$", "", a) == sub("\\.1$", "", b) & a != b
                     )))

  cors[idx]
}


plot_histograms <- function(Z, name, file) {

  cors <- get_cor_matrix(Z)

  rep_cor <- get_replicate_cor(cors)

  svg(file, width = 5, height = 6.2)

  par(mfrow = c(2, 1))

  # ---- ALL correlations ----
  h <- hist(cors,
            xlim = c(-0.4, 1),
            border = NA,
            col = "#666666",
            main = paste(name, "All correlations"),
            xlab = "Correlation")

  # ---- replicate correlations ----
  hist(rep_cor,
       xlim = c(-0.4, 1),
       breaks = h$breaks,
       col = "#386cb0",
       border = NA,
       main = paste(name, "Replicate correlations"),
       xlab = "Correlation")

  dev.off()
}

plot_histograms(Z_NTS2_NTS2_1, "NTS2", "NTS2_hists.svg")
plot_histograms(Z_NTS3_NTS3_1, "NTS3", "NTS3_hists.svg")
plot_histograms(Z_NTS4_NTS4_1, "NTS4", "NTS4_hists.svg")

## 

n_targets_list <- lapply(seq_along(upregulated), function(i) {
  data.frame(
    Sample = names(upregulated[[i]]),
    N = sapply(upregulated[[i]], length)
  )
})

## =========================
## LOAD PACKAGES
## =========================
library(topGO)
library(dplyr)
library(ggplot2)

## =========================
## CLEAN GENE FUNCTION
## =========================
clean_gene <- function(x) sub("_P.*$", "", x)

## =========================
## LOAD GO DATA
## =========================
go <- read.csv("B73_GO.out.csv", stringsAsFactors = FALSE)

go$qpid <- clean_gene(go$qpid)

go_filt <- go %>%
  filter(!is.na(ARGOT_PPV), ARGOT_PPV > 0.5)

go_filt$goid <- paste0("GO:", go_filt$goid)

gene2GO <- split(go_filt$goid, go_filt$qpid)

all_genes <- unique(go_filt$qpid)

## =========================
## MERGE REPLICATES (.1)
## =========================
clean_sample <- function(x) sub("\\.1$", "", x)

upregulated <- lapply(upregulated, function(dataset) {

  sample_names <- names(dataset)

  groups <- split(sample_names, clean_sample(sample_names))

  merged <- lapply(groups, function(samps) {
    unique(clean_gene(unlist(dataset[samps])))
  })

  merged
})

names(upregulated) <- c("NTS2", "NTS3", "NTS4")

## =========================
## TOPGO FUNCTION
## =========================
run_topGO <- function(geneList, ontology, gene2GO_list) {

  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = gene2GO_list
  )

  result <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

  GenTable(
    GOdata,
    Fishers = result,
    useLevels = TRUE,
    topNodes = 30
  )
}

## =========================
## STORAGE
## =========================
topGO_results <- list()

## =========================
## MAIN LOOP
## =========================
for (ds in names(upregulated)) {

  topGO_results[[ds]] <- list()

  for (sample in names(upregulated[[ds]])) {

    genes <- unique(upregulated[[ds]][[sample]])

    ## skip small sets
    if (length(genes) < 5) next

    ## build geneList (CRITICAL FIX)
    geneList <- factor(
      as.integer(all_genes %in% genes),
      levels = c(0, 1)
    )
    names(geneList) <- all_genes

    ## must have BOTH levels
    if (length(unique(geneList)) < 2) next

    ## run GO
    BP <- run_topGO(geneList, "BP", gene2GO)
    MF <- run_topGO(geneList, "MF", gene2GO)
    CC <- run_topGO(geneList, "CC", gene2GO)

    GO_table <- bind_rows(
      mutate(BP, Ontology = "BP"),
      mutate(MF, Ontology = "MF"),
      mutate(CC, Ontology = "CC")
    )

    GO_table$Fishers <- as.numeric(GO_table$Fishers)

    topGO_results[[ds]][[sample]] <- GO_table

    cat("Completed:", ds, "-", sample, "\n")
  }
}

## =========================
## CHECK OUTPUT
## =========================
str(topGO_results, max.level = 2)

## GO Terms across sample 
library(dplyr)

all_GO <- bind_rows(
  lapply(names(topGO_results), function(dataset){

    bind_rows(
      lapply(names(topGO_results[[dataset]]), function(sample){

        df <- topGO_results[[dataset]][[sample]]

        if(is.null(df)) return(NULL)

        df %>%
          mutate(
            Dataset = dataset,
            Sample = sample
          )

      })
    )

  })
)

save(topGO_results, all_GO, file = "GO_analysis.RData")

## keep only significant GO terms
sig_GO <- all_GO %>%
  filter(Fishers < 0.05)

library(dplyr)

## Create a unique sample identifier
all_GO <- all_GO %>%
  mutate(SampleID = paste(Dataset, Sample, sep = "_"))

## Total number of samples with GO results
samples_with_GO <- all_GO %>%
  distinct(SampleID) %>%
  nrow()

## Samples with at least one significant GO term
samples_with_sigGO <- all_GO %>%
  filter(Fishers < 0.05) %>%
  distinct(SampleID) %>%
  nrow()

## Report
cat("Samples with GO terms:", samples_with_GO, "\n") #476
cat("Samples with significant GO terms:", samples_with_sigGO, "\n") #343
cat("Total samples:", 492, "\n") #492 
cat(sprintf("GO assigned: %.1f%%\n", 100 * samples_with_GO / 492)) #96.7%
cat(sprintf("Significant GO: %.1f%%\n", 100 * samples_with_sigGO / 492)) #69.7% 

## ============================================================
## GO TERM FREQUENCY ACROSS SAMPLES
## ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

## Total number of samples
total_samples <- 492

## Create a unique sample identifier
sig_GO <- sig_GO %>%
  mutate(SampleID = paste(Dataset, Sample, sep = "_"))

## Count how many samples each GO term appears in
GO_counts <- sig_GO %>%
  group_by(GO.ID, Term, Ontology) %>%
  summarise(
    n_samples = n_distinct(SampleID),
    .groups = "drop"
  ) %>%
  mutate(
    percent_samples = round(100 * n_samples / total_samples, 1)
  ) %>%
  arrange(desc(n_samples))

## Save results
write.csv(GO_counts,
          "GO_term_frequency_across_samples.csv",
          row.names = FALSE)

## ============================================================
## SUMMARY
## ============================================================

cat("\n=============================\n")
cat("GO TERM SUMMARY\n")
cat("=============================\n")

cat("Unique significant GO terms:",
    nrow(GO_counts), "\n\n")

cat("GO terms in ALL samples:",
    sum(GO_counts$n_samples == total_samples), "\n")

cat("GO terms in >=90% samples:",
    sum(GO_counts$percent_samples >= 90), "\n")

cat("GO terms in >=75% samples:",
    sum(GO_counts$percent_samples >= 75), "\n")

cat("GO terms in >=50% samples:",
    sum(GO_counts$percent_samples >= 50), "\n\n")

cat("Top 10 most common GO terms:\n")

print(
  GO_counts %>%
    select(Term, Ontology, n_samples, percent_samples) %>%
    head(10)
)

## ============================================================
## GO TERMS FOUND IN ALL SAMPLES
## ============================================================

GO_all <- GO_counts %>%
  filter(n_samples == total_samples)

if(nrow(GO_all)==0){
  cat("\nNo GO term is significant in all", total_samples, "samples.\n")
} else{
  print(GO_all)
}

## ============================================================
## BARPLOT OF TOP 20 GO TERMS
## ============================================================

top25 <- GO_counts %>%
  slice_max(n_samples, n = 25)

p <- ggplot(top20,
            aes(reorder(Term, n_samples), n_samples)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "",
    y = "Number of samples",
    title = "Most Common Significant GO Terms"
  ) +
  theme_bw(base_size = 24)

print(p)

ggsave("Top25_GO_terms.pdf",
       p,
       width = 24,
       height = 18)

## ============================================================
## PRESENCE / ABSENCE HEATMAP
## ============================================================

## Top 50 GO terms
top_terms <- GO_counts %>%
  slice_max(n_samples, n = 50) %>%
  pull(Term)

heatmap_df <- sig_GO %>%
  filter(Term %in% top_terms) %>%
  distinct(Term, SampleID) %>%
  mutate(value = 1)

heatmap_df <- heatmap_df %>%
  complete(
    Term,
    SampleID,
    fill = list(value = 0)
  )

heatmap_plot <- ggplot(
  heatmap_df,
  aes(x = SampleID,
      y = Term,
      fill = factor(value))
) +
  geom_tile() +
  scale_fill_manual(
    values = c("0" = "white",
               "1" = "red"),
    name = "Significant"
  ) +
  labs(
    x = "Samples",
    y = "GO Term",
    title = "Presence of Significant GO Terms Across Samples"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

print(heatmap_plot)

ggsave("GO_term_heatmap.pdf",
       heatmap_plot,
       width = 32,
       height = 20)

library(dplyr)
library(tidyr)
library(ggplot2)

## ============================================================
## GO TERM HEATMAPS BY ONTOLOGY (BP / MF / CC)
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

## ============================================================
## SETTINGS
## ============================================================

top_n_terms <- 50       # number of GO terms per heatmap
ontologies <- c("BP", "MF", "CC")


## ============================================================
## MAKE HEATMAPS
## ============================================================

for (ont in ontologies) {

  cat("Creating", ont, "heatmap...\n")


  ## -----------------------------
  ## Filter ontology
  ## -----------------------------

  GO_ont <- sig_GO %>%
    filter(Ontology == ont)


  ## -----------------------------
  ## Count GO term frequency
  ## -----------------------------

  GO_counts_ont <- GO_ont %>%
    group_by(GO.ID, Term) %>%
    summarise(
      n_samples = n_distinct(SampleID),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples))


  ## -----------------------------
  ## Select top GO terms
  ## -----------------------------

  top_terms <- GO_counts_ont %>%
    slice_head(n = top_n_terms) %>%
    pull(GO.ID)

## Create presence/absence table

term_labels <- GO_ont %>%
  select(GO.ID, Term) %>%
  distinct()


heatmap_df <- GO_ont %>%
  filter(GO.ID %in% top_terms) %>%
  distinct(GO.ID, SampleID) %>%
  mutate(value = 1)


## Fill missing combinations
heatmap_df <- heatmap_df %>%
  complete(
    GO.ID,
    SampleID,
    fill = list(value = 0)
  )


## Add GO descriptions once
heatmap_df <- heatmap_df %>%
  left_join(
    term_labels,
    by = "GO.ID"
  )

  ## Order GO terms by frequency
  ## -----------------------------

  term_order <- GO_counts_ont %>%
    filter(GO.ID %in% top_terms) %>%
    arrange(n_samples) %>%
    pull(GO.ID)


  heatmap_df$GO.ID <- factor(
    heatmap_df$GO.ID,
    levels = term_order
  )


  ## -----------------------------
  ## Plot
  ## -----------------------------

  heatmap_plot <- ggplot(
    heatmap_df,
    aes(
      x = SampleID,
      y = GO.ID,
      fill = factor(value)
    )
  ) +
    geom_tile() +
    scale_fill_manual(
      values = c(
        "0" = "white",
        "1" = "red"
      ),
      name = "Significant"
    ) +
    scale_y_discrete(
      labels = function(x) {
        term_labels$Term[
          match(x, term_labels$GO.ID)
        ]
      }
    ) +
    labs(
      x = "Samples",
      y = "GO Term",
      title = paste(
        ont,
        "Significant GO Terms Across Samples"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )


  print(heatmap_plot)


  ## -----------------------------
  ## Save PDF
  ## -----------------------------

  ggsave(
    filename = paste0(
      ont,
      "_GO_term_heatmap.pdf"
    ),
    plot = heatmap_plot,
    width = 18,
    height = 12
  )


  cat("Finished", ont, "\n\n")

}

cat("All GO heatmaps completed.\n")

## ============================================================
## GO ENRICHMENT BUBBLE PLOTS BY ONTOLOGY
## ============================================================

library(dplyr)
library(ggplot2)

top_n <- 30


for (ont in c("BP", "MF", "CC")) {

  cat("Creating", ont, "bubble plot\n")


  ## Filter ontology
  GO_ont <- sig_GO %>%
    filter(Ontology == ont)


  ## Summarize GO terms
  GO_bubble <- GO_ont %>%
    group_by(GO.ID, Term) %>%
    summarise(
      n_samples = n_distinct(SampleID),
      min_Fisher = min(Fishers),
      .groups = "drop"
    ) %>%
    mutate(
      log10_p = -log10(min_Fisher),
      percent_samples = 100 * n_samples / 492
    ) %>%
    arrange(desc(n_samples)) %>%
    slice_head(n = top_n)


  ## Order by frequency using GO.ID
  GO_bubble$GO.ID <- factor(
    GO_bubble$GO.ID,
    levels = rev(GO_bubble$GO.ID)
  )


  ## Bubble plot
  p <- ggplot(
    GO_bubble,
    aes(
      x = percent_samples,
      y = GO.ID,
      size = n_samples,
      color = log10_p
    )
  ) +
    geom_point(alpha = 0.8) +
    scale_size(
      range = c(3, 12),
      name = "Samples"
    ) +
    scale_color_gradient(
      name = "-log10(Fisher)",
      low = "blue",
      high = "red"
    ) +
    scale_y_discrete(
      labels = GO_bubble$Term
    ) +
    labs(
      x = "% of samples with significant GO term",
      y = "",
      title = paste(
        ont,
        "GO Enrichment Across Samples"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
    )


  print(p)


  ggsave(
    filename = paste0(
      ont,
      "_GO_bubbleplot.pdf"
    ),
    plot = p,
    width = 10,
    height = 8
  )

############################################################################
dataset <- list(
  Z_NTS2_NTS2_1, 
  Z_NTS3_NTS3_1, 
  Z_NTS4_NTS4_1
)

## Need to average replicates (colname and colname+.1 are replicates of one another)
## Need cutoff of log10(2)for samples above that (upregulated)
## Need to figure out a way to prioritize samples to look at 
## Are there any samples with really strong effects (few genes upregulated to a high degree) and if so, how many and what does the expression look like for those 
## Are there any samples with lots of genes upregulated all together to a high degree and if so, how many and what does the expression look like for those 
## Are there any samples with little to no genes being upregulated, how many and what samples are those? 



