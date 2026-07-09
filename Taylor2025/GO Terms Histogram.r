GO Terms Histogram 
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
## Collapse replicates into a single TF ID
## (e.g. "SampleID0" and "SampleID0.1" -> "SampleID0")
## ============================================================

GO_BP_all <- GO_BP_all %>%
  mutate(TF = sub("\\.[0-9]+$", "", Sample))


## ============================================================
## Count number of significant TFs (replicate-collapsed) per GO term
## ============================================================

term_TF_counts <- GO_BP_all %>%
  distinct(Term, TF) %>%          # one row per TF-term pair, regardless of how many replicates hit it
  count(Term, sort = TRUE, name = "n_TFs")

print(head(term_TF_counts, 20))


## ============================================================
## Histogram: GO BP term vs. number of significant TFs
## ============================================================

library(ggplot2)

top_n_plot <- 50   # how many terms to show in the plot

plot_df <- term_TF_counts %>%
  slice_head(n = top_n_plot) %>%
  mutate(Term = factor(Term, levels = rev(Term)))  # so bars plot in descending order

svg("GOBP_Top100_Significant.svg", width = 25, height = 25)

ggplot(plot_df, aes(x = Term, y = n_TFs)) +
  geom_col(fill = "NavyBlue") +
  coord_flip() +
  labs(
    x = "GO Term (Biological Process)",
    y = "Number of significant TFs",
    title = "GO BP terms by number of TFs with significant enrichment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14)   # <- controls the GO term labels
  )

dev.off()

ggsave("GO_BP_term_TF_histogram.pdf", plot_df, width = 15, height =25)
