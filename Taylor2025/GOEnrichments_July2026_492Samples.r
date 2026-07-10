library(dplyr)
library(tibble)

cutoff <- log10(2)

dataset <- list(
  Z_NTS2_NTS2_1, 
  Z_NTS3_NTS3_1, 
  Z_NTS4_NTS4_1
)

Z_NTS2_NTS2_1_MC <- mc_dfs[["NTS2-NTS2"]]
Z_NTS3_NTS3_1_MC <- mc_dfs[["NTS3-NTS3"]]
Z_NTS4_NTS4_1_MC <- mc_dfs[["NTS4-NTS4"]]


############################################################
## 1. Average technical replicates
############################################################

average_replicates <- function(df){

  base_names <- sub("\\.1$", "", colnames(df))

  out <- sapply(unique(base_names), function(s){

    cols <- which(base_names == s)

    if(length(cols) == 1){
      df[, cols]
    }else{
      rowMeans(df[, cols, drop = FALSE])
    }

  })

  out <- as.data.frame(out)
  rownames(out) <- rownames(df)

  out
}

datasets_avg <- lapply(dataset, average_replicates)
names(datasets_avg) <- c("NTS2","NTS3","NTS4")

############################################################
## 2. Calculate statistics for every sample
############################################################

cutoff <- log10(2)

sample_summary <- lapply(datasets_avg, function(df){

  data.frame(

    Sample = colnames(df),

    n_up =
      colSums(df >= cutoff),

    max_expression =
      apply(df,2,max),

    mean_up_expression =
      apply(df,2,function(x){

        if(any(x >= cutoff)){
          mean(x[x >= cutoff])
        }else{
          0
        }

      }),

    median_up_expression =
      apply(df,2,function(x){

        if(any(x >= cutoff)){
          median(x[x >= cutoff])
        }else{
          0
        }

      })

  )

})

############################################################
## Combine into one table
############################################################

summary_table <- bind_rows(
  Map(function(df,name){

    df$Dataset <- name
    df

  }, sample_summary, names(sample_summary))
)

############################################################
## 3. Prioritize samples
############################################################

############################################################
## Strong effects:
## few genes, but expressed very highly
############################################################

strong_effects <-
  summary_table %>%
  arrange(desc(max_expression),
          n_up)


############################################################
## Almost no response
############################################################

little_response <-
  summary_table %>%
  arrange(n_up)

############################################################
## Show top samples
############################################################

head(strong_effects,20)
metadata_unique <- metadata %>%
  select(SampleID, Gene.model, TFomeStockID) %>%
  distinct(SampleID, .keep_all = TRUE)

nrow(summary_table)

strong_effects <- summary_table %>%
  arrange(desc(max_expression), n_up) %>%
  left_join(
    metadata_unique,
    by = c("Sample" = "SampleID")
  )

nrow(strong_effects)


head(little_response,20)

############################################################
## Counts
############################################################

cat("Samples with NO upregulated genes:", ## 12 samples 
    sum(summary_table$n_up==0),"\n")

cat("Samples with fewer than 10 upregulated genes:", ## 88 Samples 
    sum(summary_table$n_up<10),"\n")

cat("Samples with more than 100 upregulated genes:", ## 161 Samples 
    sum(summary_table$n_up>100),"\n")


cutoff <- log10(2)

datasets_upregulated <- lapply(datasets_avg, function(df) {

  # Step 1: Keep samples with at least one upregulated gene
  keep_cols <- apply(df, 2, function(x) any(x >= cutoff))
  df <- df[, keep_cols, drop = FALSE]

  # Step 2: Keep genes upregulated in at least one remaining sample
  keep_rows <- apply(df, 1, function(x) any(x >= cutoff))
  df <- df[keep_rows, , drop = FALSE]

  df
})


GeneID = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/V5GeneID.csv")
TFome = read.csv("Maize_TFome_Bulk_data _ corrected_BDN.csv")

library(dplyr)
library(tidyr)

GeneID_lookup <- GeneID %>%
  select(protein.name, all.gene.IDs) %>%
  separate_rows(all.gene.IDs, sep = " ") %>%
  rename(Gene.model = all.gene.IDs)

strong_effects <- strong_effects %>%
  left_join(
    GeneID_lookup,
    by = "Gene.model"
  )

  write.csv(strong_effects, "upregulated_genes_with_strong_effects.csv")

  ## July 8th, 2026 

############################################################
## GO TERM ENRICHMENT HEATMAP
############################################################

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

############################################################
## Load GO results
############################################################

load("GO_analysis.RData")
topGO_results 

## ============================================================
## TOP 100 BIOLOGICAL PROCESS GO TERM HEATMAP
## ============================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

load("GO_analysis.RData")
topGO_results 
## ============================================================
## Parameters
## ============================================================

p_cutoff <- 0.05
top_n <- 100


## ============================================================
## Extract BP GO enrichment results from nested list
## ============================================================

GO_BP_all <- lapply(names(topGO_results), function(plate){

  plate_results <- topGO_results[[plate]]

  lapply(names(plate_results), function(sample){

    df <- plate_results[[sample]]

    df %>%
      filter(
        Ontology == "BP",
        Fishers <= p_cutoff
      ) %>%
      mutate(
        Plate = plate,
        Sample = sample,
        negLogP = -log10(Fishers)
      )

  }) %>%
    bind_rows()

}) %>%
  bind_rows()


## Check extraction

print(head(GO_BP_all))
print(dim(GO_BP_all))


top_terms <- GO_BP_all %>%
  count(Term, sort = TRUE) %>%
  slice_head(n = top_n)


print(top_terms)


## ============================================================
## Create sample x GO term matrix
## ============================================================

heatmap_df <- GO_BP_all %>%
  filter(Term %in% top_terms$Term) %>%
  select(Sample, Term, negLogP) %>%
  group_by(Sample, Term) %>%
  summarise(
    negLogP = max(negLogP),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Term,
    values_from = negLogP,
    values_fill = 0
  )


## Convert to matrix

GO_heatmap <- heatmap_df %>%
  column_to_rownames("Sample") %>%
  as.matrix()


## ============================================================
## Order GO terms by frequency
## ============================================================

GO_heatmap <- GO_heatmap[, 
                         order(
                           colSums(GO_heatmap > 0),
                           decreasing = TRUE
                         )
]


## ============================================================
## Add sample count to GO term labels
## ============================================================

term_counts <- top_terms$n

names(term_counts) <- top_terms$Term


colnames(GO_heatmap) <- paste0(
  colnames(GO_heatmap),
  "\n(n=",
  term_counts[colnames(GO_heatmap)],
  ")"
)



## ============================================================
## Transpose matrix
## Rows = GO terms
## Columns = samples
## ============================================================

GO_heatmap <- heatmap_df %>%
  column_to_rownames("Sample") %>%
  as.matrix() %>%
  t()


## ============================================================
## Order GO terms by frequency
## ============================================================

GO_heatmap <- GO_heatmap[
  order(rowSums(GO_heatmap > 0), decreasing = TRUE),
  ]


## ============================================================
## Keep sample names as column labels
## ============================================================

colnames(GO_heatmap) <- colnames(heatmap_df)[-1]



library(viridis)

max_cap <- 4

viridis_col <- colorRamp2(
  c(0, 2, 4),
  viridis(3)
)

svg(
  "BP_GO_Top100_heatmap_Ward.D2_Clustering.svg",
  width = 50,
  height = 20
)

Heatmap(
  GO_heatmap,
  name = "-log10(Fisher P)",
  col = viridis_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90
)

dev.off()

#################### Jaccard distance 
library(vegan)
svg(
  "BP_GO_Top100_heatmap_Jaccard.svg",
  width = 50,
  height = 20
)

## Binary matrix
GO_binary <- GO_heatmap > 0

## Jaccard distances
row_dist <- vegdist(GO_binary, method = "jaccard")
col_dist <- vegdist(t(GO_binary), method = "jaccard")

## Hierarchical clustering
row_hclust <- hclust(row_dist, method = "average")
col_hclust <- hclust(col_dist, method = "average")

Heatmap(
  GO_heatmap,
  name = "-log10(Fisher P)",
  col = viridis_col,
  cluster_rows = row_hclust,
  cluster_columns = col_hclust,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90
)

dev.off()

## Important to look at dynamic range 
#row_max <- apply(GO_heatmap, 1, max)

#hist(
  #row_max,
  #breaks = 30,
  #main = "Distribution of Maximum GO Term Enrichment",
  #xlab = "Maximum -log10(Fisher P) per GO term"
#)

## Should set viridis at 4 or 5 
## Signficance Heatmap 

## TOP TERMS 
## ============================================================
## Parameters
## ============================================================

p_cutoff <- 0.05
top_n <- 100
max_cap <- 4


## ============================================================
## Extract BP GO enrichment results
## ============================================================

GO_BP_all <- lapply(names(topGO_results), function(plate){

  plate_results <- topGO_results[[plate]]

  lapply(names(plate_results), function(sample){

    df <- plate_results[[sample]]

    df %>%
      filter(
        Ontology == "BP",
        Fishers <= p_cutoff
      ) %>%
      mutate(
        Plate = plate,
        Sample = sample,
        negLogP = -log10(Fishers)
      )

  }) %>%
    bind_rows()

}) %>%
  bind_rows()


## Check

head(GO_BP_all)
dim(GO_BP_all)


## ============================================================
## Select top 100 GO terms by strongest significance
## (lowest Fisher p-value)
## ============================================================

top_terms <- GO_BP_all %>%
  group_by(Term) %>%
  summarise(
    min_Fisher = min(Fishers),
    max_negLogP = max(negLogP),
    n_samples = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(min_Fisher) %>%
  slice_head(n = top_n)


print(top_terms)


## ============================================================
## Create GO term x sample matrix
## ============================================================

heatmap_df <- GO_BP_all %>%
  filter(Term %in% top_terms$Term) %>%
  select(Sample, Term, negLogP) %>%
  group_by(Sample, Term) %>%
  summarise(
    negLogP = max(negLogP),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = negLogP,
    values_fill = 0
  )


## ============================================================
## Convert to matrix
## Rows = GO terms
## Columns = samples
## ============================================================

GO_heatmap <- heatmap_df %>%
  column_to_rownames("Term") %>%
  as.matrix()


## ============================================================
## Order GO terms by significance
## ============================================================

GO_heatmap <- GO_heatmap[
  order(rowSums(GO_heatmap > 0), decreasing = TRUE),
]


## ============================================================
## Viridis color scale capped at -log10(p)=4
## ============================================================

library(viridis)
library(circlize)

viridis_col <- colorRamp2(
  c(0, 2, max_cap),
  viridis(3)
)

svg(
  "BP_GO_Top100_By_Significance_heatmap_WardD2.svg",
  width = 50,
  height = 20
)

Heatmap( GO_heatmap, name = "-log10(Fisher P)", col = viridis_col, cluster_rows = TRUE, cluster_columns = TRUE, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2", show_row_names = TRUE, show_column_names = TRUE, column_names_rot = 90, row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 6), column_title = "Top 100 BP GO Terms by Significance", row_title = "GO Biological Processes" )
dev.off()

## ============================================================
## Heatmap
## ============================================================
library(vegan)

## ============================================================
## Jaccard clustering based on presence/absence
## ============================================================

GO_binary <- GO_heatmap > 0

row_dist <- vegdist(GO_binary, method = "jaccard")
col_dist <- vegdist(t(GO_binary), method = "jaccard")

row_hclust <- hclust(row_dist, method = "average")
col_hclust <- hclust(col_dist, method = "average")

svg(
  "BP_GO_Top100_By_Significance_heatmap_Jaccard.svg",
  width = 50,
  height = 20
)

Heatmap(
  GO_heatmap,
  name = "-log10(Fisher P)",
  col = viridis_col,

  cluster_rows = row_hclust,
  cluster_columns = col_hclust,

  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90,

  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 6),

  column_title = "Top 100 BP GO Terms by Significance",
  row_title = "GO Biological Processes"
)


dev.off()


### Moving forward with top 100 terms using ward d2 clustering ###

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)


## ============================================================
## Parameters
## ============================================================

p_cutoff <- 0.05
top_n <- 50


## ============================================================
## Extract BP GO enrichment results from nested list
## ============================================================

GO_BP_all <- lapply(names(topGO_results), function(plate){

  plate_results <- topGO_results[[plate]]

  lapply(names(plate_results), function(sample){

    df <- plate_results[[sample]]

    df %>%
      filter(
        Ontology == "BP",
        Fishers <= p_cutoff
      ) %>%
      mutate(
        Plate = plate,
        Sample = sample,
        negLogP = -log10(Fishers)
      )

  }) %>%
    bind_rows()

}) %>%
  bind_rows()


## Check extraction

print(head(GO_BP_all))
print(dim(GO_BP_all))


top_terms <- GO_BP_all %>%
  count(Term, sort = TRUE) %>%
  slice_head(n = top_n)


print(top_terms)


## ============================================================
## Create sample x GO term matrix
## ============================================================

heatmap_df <- GO_BP_all %>%
  filter(Term %in% top_terms$Term) %>%
  select(Sample, Term, negLogP) %>%
  group_by(Sample, Term) %>%
  summarise(
    negLogP = max(negLogP),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Term,
    values_from = negLogP,
    values_fill = 0
  )


## Convert to matrix

GO_heatmap <- heatmap_df %>%
  column_to_rownames("Sample") %>%
  as.matrix()


## ============================================================
## Order GO terms by frequency
## ============================================================

GO_heatmap <- GO_heatmap[, 
                         order(
                           colSums(GO_heatmap > 0),
                           decreasing = TRUE
                         )
]


## ============================================================
## Add sample count to GO term labels
## ============================================================

term_counts <- top_terms$n

names(term_counts) <- top_terms$Term


colnames(GO_heatmap) <- paste0(
  colnames(GO_heatmap),
  "\n(n=",
  term_counts[colnames(GO_heatmap)],
  ")"
)

## ============================================================
## Transpose matrix
## Rows = GO terms
## Columns = samples
## ============================================================

GO_heatmap <- heatmap_df %>%
  column_to_rownames("Sample") %>%
  as.matrix() %>%
  t()


## ============================================================
## Order GO terms by frequency
## ============================================================

GO_heatmap <- GO_heatmap[
  order(rowSums(GO_heatmap > 0), decreasing = TRUE),
  ]


## ============================================================
## Keep sample names as column labels
## ============================================================

colnames(GO_heatmap) <- colnames(heatmap_df)[-1]

max_cap <- 4

viridis_col <- colorRamp2(
  c(0, 2, 4),
  viridis(3)
)

svg(
  "BP_GO_Top50_heatmap_Ward.D2_Clustering.svg",
  width = 50,
  height = 20
)

Heatmap(
  GO_heatmap,
  name = "-log10(Fisher P)",
  col = viridis_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90
)

dev.off()

## Looking at Flowwering time set of TFs 

S034674-02 
S034674-163
S036682-65 
S036682-125 
S036682-27 
S034674-55
S036682-50
S034674-56 
S034674-232
S036682-156
S036682-228
S034674-22
S036682-174
S031591-148 
S034674-21 
S034674-269
S034674-233
S034674-04
S036682-220
S036682-25
S034674-218
S036682-144
S036682-19
S036682-94
S031591-134
S036682-199

## SampleIDs of interest
sample_ids <- c(
  "S034674-02",
  "S034674-163",
  "S036682-65",
  "S036682-125",
  "S036682-27",
  "S034674-55",
  "S036682-50",
  "S034674-56",
  "S034674-232",
  "S036682-156",
  "S036682-228",
  "S034674-22",
  "S036682-174",
  "S031591-148",
  "S034674-21",
  "S034674-269",
  "S034674-233",
  "S034674-04",
  "S036682-220",
  "S036682-25",
  "S034674-218",
  "S036682-144",
  "S036682-19",
  "S036682-94",
  "S031591-134",
  "S036682-199"
)

## Match metadata
sample_metadata <- metadata %>%
  dplyr::filter(SampleID %in% sample_ids) %>%
  dplyr::select(SampleID, Gene.model, TFomeStockID)

## Check for missing IDs
missing_ids <- setdiff(sample_ids, sample_metadata$SampleID)
missing_ids

## Expand all.gene.IDs so each gene ID gets its own row
GeneID_lookup <- GeneID %>%
  dplyr::select(protein.name, all.gene.IDs) %>%
  tidyr::separate_rows(all.gene.IDs, sep = "\\s+") %>%
  dplyr::rename(Gene.model = all.gene.IDs) %>%
  dplyr::distinct()

## Join to your sample metadata
sample_metadata <- sample_metadata %>%
  left_join(GeneID_lookup, by = "Gene.model")

## Write CSV
write.csv(
  sample_metadata,
  "Selected_SampleIDs_metadata.csv",
  row.names = FALSE
)

library(dplyr)
library(ggplot2)

sample_metadata <- sample_metadata %>%
  mutate(
    TF_family = protein.name %>%
      sub("^Zm", "", .) %>%      # Remove "Zm" prefix
      sub("[0-9].*$", "", .)     # Remove numbers and everything after
  )

## Count families
family_counts <- sample_metadata %>%
  count(TF_family, sort = TRUE)


library(ggplot2)

svg(
  "TFDistributioninFloweringCluster.svg",
  width = 25,
  height = 20
)

ggplot(family_counts,
       aes(x = reorder(TF_family, n), y = n)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, max(family_counts$n), by = 1)) +
  labs(x = "TF Group", y = "Number of TFs in Flowering Cluster") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 22),
    axis.ticks = element_line(linewidth = 0.8)
  )

dev.off()