## ============================================================
## Tissue Specific Expression Heatmap on the Samples in the
## Flowering Time Cluster — labeled with protein.name
## ============================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments")

## ---- 1. Load gene ID -> protein name reference table ---------------------
## (must be loaded BEFORE it's used below)
GeneID <- read.csv("V5GeneID.csv", check.names = FALSE, stringsAsFactors = FALSE)

## ---- 2. Expand all.gene.IDs into one row per individual ID -----------------
gene_id_long <- GeneID %>%
  select(protein.name, all.gene.IDs) %>%
  separate_rows(all.gene.IDs, sep = "\\s+") %>%   # split on whitespace
  filter(all.gene.IDs != "")                       # drop empty strings

## ---- 3. Keep only the v4-style IDs (Zm00001d...) to match your expression data ----
gene_id_v4 <- gene_id_long %>%
  filter(str_starts(all.gene.IDs, "Zm00001d"))

## ---- 4. Check for duplicates before building the lookup ---------------------
dup_check <- gene_id_v4 %>% count(all.gene.IDs) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  message(nrow(dup_check), " v4 ID(s) map to more than one protein.name — inspect these:")
  print(gene_id_v4 %>% filter(all.gene.IDs %in% dup_check$all.gene.IDs))
  ## keep first match only, so the lookup stays 1:1
  gene_id_v4 <- gene_id_v4 %>% distinct(all.gene.IDs, .keep_all = TRUE)
}

## ---- 5. Build the lookup: v4 gene ID -> protein name -------------------------
id_to_protein <- setNames(gene_id_v4$protein.name, gene_id_v4$all.gene.IDs)
head(id_to_protein)

## ---- 6. Apply to your transposed matrix (columns = genes) ----------------------
## mat_scaled_t must already exist from your earlier heatmap-building steps,
## with gene IDs (Zm00001d...) as column names
current_gene_ids <- colnames(mat_scaled_t)

## sanity check the match rate before renaming
message(sum(current_gene_ids %in% names(id_to_protein)), " of ",
        length(current_gene_ids), " gene IDs matched to a protein name.")
unmatched <- current_gene_ids[!(current_gene_ids %in% names(id_to_protein))]
if (length(unmatched) > 0) {
  message("Note: kept original ID for unmatched gene(s): ",
          paste(unmatched, collapse = ", "))
}

new_col_names <- ifelse(
  current_gene_ids %in% names(id_to_protein),
  id_to_protein[current_gene_ids],
  current_gene_ids   # fall back to original ID if no match found
)

colnames(mat_scaled_t) <- new_col_names

## ---- 7. Re-draw heatmap with new column labels ----------------------------------
ht <- Heatmap(
  mat_scaled_t,
  name = "z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 90,
  column_title = "Flowering Time Genes: Tissue Expression (row z-score)"
)

svg("maize_floweringTime_heatmap.svg", width = 8, height = 24)
draw(ht)
dev.off()
