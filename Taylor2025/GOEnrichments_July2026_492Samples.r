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


  ..$ S036682-220:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:2000028" "GO:0048573" "GO:0010228" "GO:0009648" ...
  .. ..$ Term       : chr [1:90] "regulation of photoperiodism, flowering" "photoperiodism, flowering" "vegetative to reproductive phase transit..." "photoperiodism" ...
  .. ..$ Level      : int [1:90] 10 9 8 6 4 7 6 4 5 5 ...
  .. ..$ Annotated  : int [1:90] 7 7 7 7 17 17 17 21 23 23 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0.01 0.01 0.01 0.01 0.03 0.03 0.03 0.04 0.05 0.05 ...
  .. ..$ Fishers    : num [1:90] 0.014 0.014 0.014 0.014 0.033 0.033 0.033 0.041 0.045 0.045 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-221:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902356" "GO:1901001" "GO:0048585" "GO:0009651" ...
  .. ..$ Term       : chr [1:90] "oxaloacetate(2-) transmembrane transport" "negative regulation of response to salt ..." "negative regulation of response to stimu..." "response to salt stress" ...
  .. ..$ Level      : int [1:90] 11 8 5 5 6 7 4 4 3 5 ...
  .. ..$ Annotated  : int [1:90] 1 6 25 50 50 50 50 63 98 145 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0 0.01 0.05 0.1 0.1 0.1 0.1 0.12 0.19 0.29 ...
  .. ..$ Fishers    : num [1:90] 0.002 0.012 0.049 0.096 0.096 0.096 0.096 0.12 0.184 0.266 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-222:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902395" "GO:1902348" "GO:1901600" "GO:0006714" ...
  .. ..$ Term       : chr [1:90] "regulation of 1-deoxy-D-xylulose-5-phosp..." "cellular response to strigolactone" "strigolactone metabolic process" "sesquiterpenoid metabolic process" ...
  .. ..$ Level      : int [1:90] 6 7 9 8 5 6 7 6 5 4 ...
  .. ..$ Annotated  : int [1:90] 2 2 7 7 8 8 10 11 18 22 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0 0 0.01 0.01 0.02 0.02 0.02 0.02 0.04 0.04 ...
  .. ..$ Fishers    : num [1:90] 0.0039 0.0039 0.0138 0.0138 0.0157 0.0157 0.0196 0.0216 0.0352 0.0429 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-223:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:0010150" "GO:0090693" "GO:1900055" "GO:0072593" ...
  .. ..$ Term       : chr [1:90] "leaf senescence" "plant organ senescence" "regulation of leaf senescence" "reactive oxygen species metabolic proces..." ...
  .. ..$ Level      : int [1:90] 9 7 10 4 6 8 7 6 5 8 ...
  .. ..$ Annotated  : int [1:90] 7 7 7 9 9 9 9 10 17 21 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0.02 0.02 0.02 0.03 0.03 0.03 0.03 0.03 0.05 0.06 ...
  .. ..$ Fishers    : num [1:90] 0.021 0.021 0.021 0.026 0.026 0.026 0.026 0.029 0.05 0.061 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-224:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902356" "GO:1903602" "GO:0009310" "GO:0006598" ...
  .. ..$ Term       : chr [1:90] "oxaloacetate(2-) transmembrane transport" "thermospermine catabolic process" "amine catabolic process" "polyamine catabolic process" ...
  .. ..$ Level      : int [1:90] 11 9 6 8 7 8 7 6 6 4 ...
  .. ..$ Annotated  : int [1:90] 1 4 4 4 4 5 5 7 64 64 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 2 2 ...
  .. ..$ Expected   : num [1:90] 0 0.02 0.02 0.02 0.02 0.02 0.02 0.03 0.32 0.32 ...
  .. ..$ Fishers    : num [1:90] 0.0049 0.0196 0.0196 0.0196 0.0196 0.0245 0.0245 0.0341 0.0347 0.0347 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-225:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902356" "GO:1902395" "GO:1902600" "GO:0098662" ...
  .. ..$ Term       : chr [1:90] "oxaloacetate(2-) transmembrane transport" "regulation of 1-deoxy-D-xylulose-5-phosp..." "proton transmembrane transport" "inorganic cation transmembrane transport" ...
  .. ..$ Level      : int [1:90] 11 6 8 7 7 6 6 5 6 5 ...
  .. ..$ Annotated  : int [1:90] 1 2 109 126 126 127 129 130 150 274 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 2 ...
  .. ..$ Expected   : num [1:90] 0 0.01 0.32 0.37 0.37 0.38 0.38 0.38 0.44 0.81 ...
  .. ..$ Fishers    : num [1:90] 0.003 0.0059 0.2893 0.3287 0.3287 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-226:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1904430" "GO:0002832" "GO:0032102" "GO:0031348" ...
  .. ..$ Term       : chr [1:90] "negative regulation of t-circle formatio..." "negative regulation of response to bioti..." "negative regulation of response to exter..." "negative regulation of defense response" ...
  .. ..$ Level      : int [1:90] 12 6 6 7 9 9 10 7 9 8 ...
  .. ..$ Annotated  : int [1:90] 1 3 3 3 3 4 4 5 5 5 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0.01 0.02 0.02 0.02 0.02 0.03 0.03 0.04 0.04 0.04 ...
  .. ..$ Fishers    : num [1:90] 0.0079 0.0235 0.0235 0.0235 0.0235 0.0312 0.0312 0.0389 0.0389 0.0389 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-227:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902356" "GO:0051049" "GO:0032879" "GO:1905952" ...
  .. ..$ Term       : chr [1:90] "oxaloacetate(2-) transmembrane transport" "regulation of transport" "regulation of localization" "regulation of lipid localization" ...
  .. ..$ Level      : int [1:90] 11 5 4 5 6 8 9 8 9 7 ...
  .. ..$ Annotated  : int [1:90] 1 20 20 2 2 2 2 2 2 2 ...
  .. ..$ Significant: int [1:90] 1 2 2 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0.01 0.18 0.18 0.02 0.02 0.02 0.02 0.02 0.02 0.02 ...
  .. ..$ Fishers    : num [1:90] 0.0089 0.0123 0.0123 0.0177 0.0177 0.0177 0.0177 0.0177 0.0177 0.0177 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-228:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1902356" "GO:0048573" "GO:2000028" "GO:0009648" ...
  .. ..$ Term       : chr [1:90] "oxaloacetate(2-) transmembrane transport" "photoperiodism, flowering" "regulation of photoperiodism, flowering" "photoperiodism" ...
  .. ..$ Level      : int [1:90] 11 9 10 6 8 9 10 8 3 4 ...
  .. ..$ Annotated  : int [1:90] 1 7 7 7 7 8 8 12 98 228 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 2 3 ...
  .. ..$ Expected   : num [1:90] 0 0.03 0.03 0.03 0.03 0.04 0.04 0.06 0.48 1.12 ...
  .. ..$ Fishers    : num [1:90] 0.0049 0.0341 0.0341 0.0341 0.0341 0.0389 0.0389 0.0579 0.0762 0.0783 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  ..$ S036682-229:'data.frame': 90 obs. of  8 variables:
  .. ..$ GO.ID      : chr [1:90] "GO:1903602" "GO:0042402" "GO:0009310" "GO:0006598" ...
  .. ..$ Term       : chr [1:90] "thermospermine catabolic process" "cellular biogenic amine catabolic proces..." "amine catabolic process" "polyamine catabolic process" ...
  .. ..$ Level      : int [1:90] 9 7 6 8 8 7 6 7 5 6 ...
  .. ..$ Annotated  : int [1:90] 4 4 4 4 5 5 7 8 8 8 ...
  .. ..$ Significant: int [1:90] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Expected   : num [1:90] 0.01 0.01 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 ...
  .. ..$ Fishers    : num [1:90] 0.012 0.012 0.012 0.012 0.015 0.015 0.021 0.024 0.024 0.024 ...
  .. ..$ Ontology   : chr [1:90] "BP" "BP" "BP" "BP" ...
  .. [list output truncated]
NULL
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


### 
### Moving forward with broad GO terms (Level <= 5) using Ward D2 clustering ###

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

## GO hierarchy level cutoff
## Lower = broader categories
GO_level_cutoff <- 5


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


## ============================================================
## Collapse to broader GO categories
## ============================================================

GO_BP_all <- GO_BP_all %>%
  filter(
    Level <= GO_level_cutoff
  )


print(dim(GO_BP_all))


## ============================================================
## Select top 50 broad GO terms by sample frequency
## ============================================================

top_terms <- GO_BP_all %>%
  count(Term, sort = TRUE) %>%
  slice_head(n = top_n)


print(top_terms)



## ============================================================
## Create sample x GO term matrix
## ============================================================

heatmap_df <- GO_BP_all %>%
  filter(
    Term %in% top_terms$Term
  ) %>%
  select(
    Sample,
    Term,
    negLogP
  ) %>%
  group_by(
    Sample,
    Term
  ) %>%
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
## Order GO terms by frequency
## ============================================================

GO_heatmap <- GO_heatmap[
  order(
    rowSums(GO_heatmap > 0),
    decreasing = TRUE
  ),
]



## ============================================================
## Color scale capped at -log10(P)=4
## ============================================================

max_cap <- 4

viridis_col <- colorRamp2(
  c(0, 2, 4),
  viridis(3)
)



## ============================================================
## Ward.D2 clustered heatmap
## ============================================================

svg(
  "BP_GO_Broad_Level5_Top50_Ward.D2.svg",
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

  column_names_rot = 90,

  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 6),

  column_title = paste0(
    "Top ",
    top_n,
    " Broad BP GO Terms (GO Level ≤ ",
    GO_level_cutoff,
    ")"
  ),

  row_title = "Biological Processes"
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