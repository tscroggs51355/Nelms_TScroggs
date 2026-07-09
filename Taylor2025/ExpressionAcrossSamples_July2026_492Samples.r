setwd("C:/Users/taylo/Desktop/TF Project Master Doucments")

library(dplyr)


## ============================================================
## Load metadata
## ============================================================

metadata <- read.csv(
    "Metadata_FullPlate_Jan2026.csv",
    check.names = FALSE
)



## ============================================================
## Load normalized expression files
## ============================================================

Z_NTS2_NTS2_1 <- read.csv(
    "Z_NTS2_NTS2_1_Normalized_Unfiltered.csv",
    row.names = 1,
    check.names = FALSE
)

Z_NTS3_NTS3_1 <- read.csv(
    "Z_NTS3_NTS3_1_Normalized_Unfiltered.csv",
    row.names = 1,
    check.names = FALSE
)

Z_NTS4_NTS4_1 <- read.csv(
    "Z_NTS4_NTS4_1_Normalized_Unfiltered.csv",
    row.names = 1,
    check.names = FALSE
)



## ============================================================
## Rename wells to SampleID using metadata
## ============================================================

rename_to_sampleID <- function(df, plate_name){

    meta_sub <- metadata %>%
        filter(Assay.Plate == plate_name)


    lookup <- setNames(
        meta_sub$SampleID,
        meta_sub$Array.Well
    )


    old_names <- colnames(df)


    new_names <- ifelse(
        old_names %in% names(lookup),
        lookup[old_names],
        old_names
    )


    colnames(df) <- new_names

    df
}



Z_NTS2_NTS2_1 <- rename_to_sampleID(
    Z_NTS2_NTS2_1,
    "NTS2-NTS2"
)


Z_NTS3_NTS3_1 <- rename_to_sampleID(
    Z_NTS3_NTS3_1,
    "NTS3-NTS3"
)


Z_NTS4_NTS4_1 <- rename_to_sampleID(
    Z_NTS4_NTS4_1,
    "NTS4-NTS4"
)



## ============================================================
## Rename every mCherry column individually
## ============================================================

fix_mCherry_names <- function(df, plate){

    mc <- which(colnames(df) == "mCherry")


    if(length(mc) > 0){

        colnames(df)[mc] <- paste0(
            plate,
            "_MC_",
            seq_along(mc)
        )

    }

    df
}



Z_NTS2_NTS2_1 <- fix_mCherry_names(
    Z_NTS2_NTS2_1,
    "NTS2"
)


Z_NTS3_NTS3_1 <- fix_mCherry_names(
    Z_NTS3_NTS3_1,
    "NTS3"
)


Z_NTS4_NTS4_1 <- fix_mCherry_names(
    Z_NTS4_NTS4_1,
    "NTS4"
)



## ============================================================
## Split mCherry from experimental samples
## ============================================================

separate_mCherry <- function(df){

    mc <- grep(
        "_MC_",
        colnames(df),
        value = TRUE
    )


    list(

        samples =
            df[, !colnames(df) %in% mc, drop = FALSE],

        mCherry =
            df[, mc, drop = FALSE]

    )
}



NTS2_split <- separate_mCherry(Z_NTS2_NTS2_1)

NTS3_split <- separate_mCherry(Z_NTS3_NTS3_1)

NTS4_split <- separate_mCherry(Z_NTS4_NTS4_1)



## ============================================================
## Average only SampleID technical replicates
## SampleID and SampleID.1
## ============================================================

average_replicates <- function(df){


    base <- sub(
        "\\.1$",
        "",
        colnames(df)
    )


    groups <- split(
        seq_along(base),
        base
    )


    averaged <- sapply(
        groups,
        function(x){

            rowMeans(
                df[, x, drop = FALSE],
                na.rm = TRUE
            )

        }
    )


    averaged <- as.data.frame(averaged)


    colnames(averaged) <- names(groups)


    averaged

}



Z_NTS2_avg <- average_replicates(
    NTS2_split$samples
)


Z_NTS3_avg <- average_replicates(
    NTS3_split$samples
)


Z_NTS4_avg <- average_replicates(
    NTS4_split$samples
)



## ============================================================
## Combine mCherry controls WITHOUT averaging
## ============================================================

NTS2_MC <- NTS2_split$mCherry

NTS3_MC <- NTS3_split$mCherry

NTS4_MC <- NTS4_split$mCherry



mCherry_controls <- cbind(
    NTS2_MC,
    NTS3_MC,
    NTS4_MC
)



## ============================================================
## Find common genes
## ============================================================

common_genes <- Reduce(
    intersect,
    list(
        rownames(Z_NTS2_avg),
        rownames(Z_NTS3_avg),
        rownames(Z_NTS4_avg),
        rownames(mCherry_controls)
    )
)



## ============================================================
## Final expression objects
## ============================================================


Z_all_avg <- cbind(

    Z_NTS2_avg[common_genes, , drop = FALSE],

    Z_NTS3_avg[common_genes, , drop = FALSE],

    Z_NTS4_avg[common_genes, , drop = FALSE]

)



mCherry_controls <- mCherry_controls[
    common_genes,
    ,
    drop = FALSE
]



## ============================================================
## Checks
## ============================================================


print(dim(Z_all_avg))

print(dim(mCherry_controls))


print(head(colnames(Z_all_avg)))

print(head(colnames(mCherry_controls)))



## Check mCherry are unique

print(
    cor(
        mCherry_controls[,1],
        mCherry_controls[,2],
        use = "complete.obs"
    )
)

## Top 100 Most Variable 
## ============================================================
## Top 100 variable genes standard clustered heatmap
## ============================================================

library(ComplexHeatmap)
library(circlize)


## ------------------------------------------------------------
## Combine samples + mCherry controls
## ------------------------------------------------------------

expr_samples <- Z_all_avg[, sample_ids, drop = FALSE]


common_genes <- intersect(
    rownames(expr_samples),
    rownames(mCherry_controls)
)


expr_heatmap <- cbind(
    expr_samples[common_genes, , drop = FALSE],
    mCherry_controls[common_genes, , drop = FALSE]
)



## ------------------------------------------------------------
## Select top 100 most variable genes
## ------------------------------------------------------------

gene_var <- apply(
    expr_heatmap,
    1,
    var,
    na.rm = TRUE
)


top100 <- names(
    sort(
        gene_var,
        decreasing = TRUE
    )[1:100]
)


expr_top100 <- expr_heatmap[
    top100,
    ,
    drop = FALSE
]



## ------------------------------------------------------------
## Heatmap colors
## ------------------------------------------------------------

hmcols2 <- colorRamp2(
    c(
        min(expr_top100, na.rm = TRUE),
        0,
        max(expr_top100, na.rm = TRUE)
    ),
    c(
        "blue",
        "white",
        "red"
    )
)



## ------------------------------------------------------------
## Standard hierarchical clustering heatmap
## ------------------------------------------------------------
svg(
  "Expression_Flowering_Top100VariableGenes.svg",
  width = 50,
  height = 25
)
Heatmap(
    expr_top100,

    name = "Expression",

    col = hmcols2,

    cluster_rows = TRUE,

    cluster_columns = TRUE,

    clustering_distance_rows = "euclidean",

    clustering_distance_columns = "euclidean",

    clustering_method_rows = "complete",

    clustering_method_columns = "complete",

    show_row_names = TRUE,

    show_column_names = TRUE,

    row_names_gp = grid::gpar(
        fontsize = 6
    ),

    column_names_gp = grid::gpar(
        fontsize = 8
    ),

    column_title = "Top 100 Variable Genes"
)
dev.off()

## Top 100 most induced genees 
## Calculate maximum induction across samples

gene_max <- apply(
    expr_samples,
    1,
    max,
    na.rm = TRUE
)


## Select top 100 most induced genes

top100_upregulated <- names(
    sort(
        gene_max,
        decreasing = TRUE
    )[1:100]
)


expr_top100 <- expr_heatmap[
    top100_upregulated,
    ,
    drop = FALSE
]

## ============================================================
## Replace SampleID column names with protein.name
## ============================================================

protein_lookup <- sample_metadata %>%
    dplyr::select(SampleID, protein.name) %>%
    dplyr::distinct()


new_names <- protein_lookup$protein.name[
    match(
        colnames(expr_top100),
        protein_lookup$SampleID
    )
]


## Keep original names if no match (mCherry controls, etc.)
new_names[is.na(new_names)] <- colnames(expr_top100)[is.na(new_names)]


colnames(expr_top100) <- new_names

## ------------------------------------------------------------
## Standard hierarchical clustering heatmap
## ------------------------------------------------------------
svg(
  "Expression_Flowering_Top100_BasedonExpression_Scaling_Clustering_Columns.svg",
  width = 40,
  height = 25
)

hmcols2 <- colorRamp2(
    c(
        min(expr_top100, na.rm = TRUE),
        0,
        max(expr_top100, na.rm = TRUE)
    ),
    c(
        "blue",
        "white",
        "red"
    )
)

Heatmap(
    expr_top100,

    name = "Expression",

    col = hmcols2,

    cluster_rows = TRUE,

    cluster_columns = TRUE,

    clustering_distance_rows = "euclidean",

    clustering_distance_columns = "euclidean",

    clustering_method_rows = "complete",

    clustering_method_columns = "complete",

    show_row_names = TRUE,

    show_column_names = TRUE,

    row_names_gp = grid::gpar(
        fontsize = 6
    ),

    column_names_gp = grid::gpar(
        fontsize = 8
    ),

    column_title = "Top 100 Most Highly Expresed Genes"
)
dev.off()

