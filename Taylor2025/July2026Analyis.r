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

## keep only significant GO terms
sig_GO <- all_GO %>%
  filter(Fishers < 0.05)

GO_frequency <- sig_GO %>%
  group_by(GO.ID, Term, Ontology) %>%
  summarise(
    n_samples = n(),
    .groups="drop"
  ) %>%
  arrange(desc(n_samples))

  library(ggplot2)

GO_frequency %>%
  slice_max(n_samples, n=20) %>%
  ggplot(aes(reorder(Term, n_samples), n_samples))+
  geom_col(fill="steelblue")+
  coord_flip()+
  labs(
    x="GO Term",
    y="Number of Samples",
    title="Most Frequently Enriched GO Terms"
  )+
  theme_bw()

  ## Unique Go Function 
  unique_GO <- sig_GO %>%
  group_by(GO.ID, Term, Dataset) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(GO.ID, Term) %>%
  mutate(
    datasets_present=n()
  ) %>%
  filter(datasets_present==1)
  
  unique_GO %>%
  arrange(Dataset, desc(n))

  ## Heatmap of presence/absence 
  heatmap_df <- sig_GO %>%
  mutate(value=1) %>%
  select(Term, Sample, value)

ggplot(heatmap_df,
       aes(Sample, Term, fill=value))+
    geom_tile()+
    scale_fill_gradient(low="white", high="darkred")+
    theme_bw()

## Counting GO terms 
GO_counts <- sig_GO %>%
  mutate(SampleID = paste(Dataset, Sample, sep = "_")) %>%
  group_by(GO.ID, Term, Ontology) %>%
  summarise(
    n_samples = n_distinct(SampleID),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

  head(GO_counts, 20)

  ## Plotting distribution across samples 
  library(ggplot2)

ggplot(GO_counts, aes(x = n_samples)) +
  geom_histogram(binwidth = 1, color = "black", fill = "steelblue") +
  scale_x_continuous(breaks = seq(1, max(GO_counts$n_samples), by = 1)) +
  labs(
    x = "Number of samples containing the GO term",
    y = "Number of GO terms",
    title = "Distribution of GO terms across samples"
  ) +
  theme_bw()

GO_counts %>%
  filter(n_samples >= 10) %>%
  arrange(desc(n_samples))

  total_samples <- 492

GO_summary <- sig_GO %>%
  mutate(SampleID = paste(Dataset, Sample, sep = "_")) %>%   # unique sample ID
  group_by(GO.ID, Term, Ontology) %>%
  summarise(
    n_samples = n_distinct(SampleID),
    percent_samples = round(100 * n_samples / total_samples, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))