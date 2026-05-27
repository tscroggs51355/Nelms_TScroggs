#######################################################################################3

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments") 

metadata = read.csv("Metadata_FullPlate_Jan2026.csv")
### 1. Extract quadrant from SampleName
metadata$Quadrant <- sub(".*-([A-D])_.*", "\\1", metadata$SampleName)

### 2. Build metadata lookup
meta_lookup <- metadata[, c("Assay.Plate", "Array.Well", "SampleID", "SampleName", "Quadrant")]

A2_NTS2_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS2_NTS2_1.csv", row.names = 1, check.names = FALSE)
A2_NTS3_NTS3_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS3_NTS3_1.csv", row.names = 1, check.names = FALSE)
A2_NTS4_NTS4_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS4_NTS4_1_PreNormalization.csv", row.names = 1, check.names = FALSE)
A2_NTS1_NTS2_1 = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/PreNormalization/A2_NTS1_NTS2_1_PreNormalization.csv", row.names = 1, check.names = FALSE)

controls_by_plate <- split(
  metadata$Array.Well[metadata$SampleID == "mCherry"],
  metadata$Assay.Plate[metadata$SampleID == "mCherry"]
)

Z_NTS2_NTS2_1 = read.csv("Z_NTS2_NTS2_1_Normalized_Unfiltered.csv")
Z_NTS3_NTS3_1 = read.csv("Z_NTS3_NTS3_1_Normalized_Unfiltered.csv")
Z_NTS4_NTS4_1 = read.csv("Z_NTS4_NTS4_1_Normalized_Unfiltered.csv")
Z_NTS1_NTS2_1 = read.csv("Z_NTS1_NTS2_1_Normalized_Unfiltered.csv")


Z_NTS2_NTS2_1 <- dfs$Z_NTS2_NTS2_1
Z_NTS3_NTS3_1 <- dfs$Z_NTS3_NTS3_1
Z_NTS4_NTS4_1 <- dfs$Z_NTS4_NTS4_1
Z_NTS1_NTS2_1 <- dfs$Z_NTS1_NTS2_1



Z_NTS2_NTS2_1_filtered <- Z_NTS2_NTS2_1[rowSums(abs(Z_NTS2_NTS2_1) >= log10(2)) > 0, ]
Z_NTS3_NTS3_1_filtered <- Z_NTS3_NTS3_1[rowSums(abs(Z_NTS3_NTS3_1) >= log10(2)) > 0, ]
Z_NTS4_NTS4_1_filtered <- Z_NTS4_NTS4_1[rowSums(abs(Z_NTS4_NTS4_1) >= log10(2)) > 0, ]
Z_NTS1_NTS2_1_filtered <- Z_NTS1_NTS2_1[rowSums(abs(Z_NTS1_NTS2_1) >= log10(2)) > 0, ]

write.csv(Z_NTS2_NTS2_1_filtered, "Z_NTS2_NTS2_1_filtered.csv")
write.csv(Z_NTS3_NTS3_1_filtered, "Z_NTS3_NTS3_1_filtered.csv")
write.csv(Z_NTS4_NTS4_1_filtered, "Z_NTS4_NTS4_1_filtered.csv")
write.csv(Z_NTS1_NTS2_1_filtered, "Z_NTS1_NTS2_1_filtered.csv")



dfs <- list(
  Z_NTS2_NTS2_1 = Z_NTS2_NTS2_1,
  Z_NTS3_NTS3_1 = Z_NTS3_NTS3_1,
  Z_NTS4_NTS4_1 = Z_NTS4_NTS4_1,
  Z_NTS1_NTS2_1 = Z_NTS1_NTS2_1,
  A2_NTS2_NTS2_1 = A2_NTS2_NTS2_1,
  A2_NTS3_NTS3_1 = A2_NTS3_NTS3_1,
  A2_NTS4_NTS4_1 = A2_NTS4_NTS4_1,
  A2_NTS1_NTS2_1 = A2_NTS1_NTS2_1, 
  Z_NTS2_NTS2_1_filtered = Z_NTS2_NTS2_1_filtered, 
  Z_NTS3_NTS3_1_filtered = Z_NTS3_NTS3_1_filtered, 
  Z_NTS4_NTS4_1_filtered = Z_NTS4_NTS4_1_filtered, 
  Z_NTS1_NTS2_1_filtered = Z_NTS1_NTS2_1_filtered
)

meta_split <- split(meta_lookup, meta_lookup$Assay.Plate)
rename_plate <- function(df, name, meta_split) {

plate <- name |>
  sub("^[^_]+_", "", x = _) |>
  sub("_filtered$", "", x = _) |>
  sub("_1$", "", x = _) |>
  gsub("_", "-", x = _)

  meta <- meta_split[[plate]]

  wells <- colnames(df)

  idx <- match(wells, meta$Array.Well)
  new_names <- meta$SampleID[idx]

  # fallback if no match
  new_names[is.na(new_names)] <- wells[is.na(new_names)]

  colnames(df) <- new_names

  df
}

dfs_renamed <- lapply(names(dfs), function(nm) {
  rename_plate(dfs[[nm]], nm, meta_split)
})
names(dfs_renamed) <- names(dfs)

Z_NTS2_NTS2_1_SID <- dfs_renamed$Z_NTS2_NTS2_1
Z_NTS3_NTS3_1_SID <- dfs_renamed$Z_NTS3_NTS3_1
Z_NTS4_NTS4_1_SID <- dfs_renamed$Z_NTS4_NTS4_1
Z_NTS1_NTS2_1_SID <- dfs_renamed$Z_NTS1_NTS2_1


A2_NTS2_NTS2_1_SID <- dfs_renamed$A2_NTS2_NTS2_1
A2_NTS3_NTS3_1_SID <- dfs_renamed$A2_NTS3_NTS3_1
A2_NTS4_NTS4_1_SID <- dfs_renamed$A2_NTS4_NTS4_1
A2_NTS1_NTS2_1_SID <- dfs_renamed$A2_NTS1_NTS2_1


  Z_NTS2_NTS2_1_filtered_SID = dfs_renamed$Z_NTS2_NTS2_1_filtered
  Z_NTS3_NTS3_1_filtered_SID = dfs_renamed$Z_NTS3_NTS3_1_filtered 
  Z_NTS4_NTS4_1_filtered_SID = dfs_renamed$Z_NTS4_NTS4_1_filtered 
  Z_NTS1_NTS2_1_filtered_SID = dfs_renamed$Z_NTS1_NTS2_1_filtered


### UMIS per Quadrant Per Plate 

datasets_A2 <- list(
  "NTS2-NTS2" = A2_NTS2_NTS2_1,
  "NTS3-NTS3" = A2_NTS3_NTS3_1,
  "NTS4-NTS4" = A2_NTS4_NTS4_1,
  "NTS1-NTS2" = A2_NTS1_NTS2_1
)

report_A2 <- do.call(rbind, lapply(names(datasets_A2), function(plate_name) {

  A2 <- datasets_A2[[plate_name]]

  meta_sub <- meta_lookup[meta_lookup$Assay.Plate == plate_name & !is.na(meta_lookup$SampleID), ]

  # safer quadrant extraction
  meta_sub$Quadrant <- regmatches(
    meta_sub$SampleName,
    regexpr("[ABCD]", meta_sub$SampleName)
  )

  quad_means <- sapply(c("A", "B", "C", "D"), function(q) {

    wells <- meta_sub$Array.Well[meta_sub$Quadrant == q]

    wells <- intersect(wells, colnames(A2))  # safer than %in%

    if (length(wells) == 0) return(NA)

    mean(colSums(A2[, wells, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
  })

  data.frame(
    Dataset    = plate_name,
    Quadrant_A = round(quad_means["A"], 2),
    Quadrant_B = round(quad_means["B"], 2),
    Quadrant_C = round(quad_means["C"], 2),
    Quadrant_D = round(quad_means["D"], 2),
    row.names  = NULL
  )
}))

write.csv(report_A2, "quadrant_avg_UMI_report_V2.csv", row.names = FALSE)

### Genes Detected per Sample 
datasets_A2 <- list(
  "NTS2-NTS2" = A2_NTS2_NTS2_1,
  "NTS3-NTS3" = A2_NTS3_NTS3_1,
  "NTS4-NTS4" = A2_NTS4_NTS4_1,
  "NTS1-NTS2" = A2_NTS1_NTS2_1
)

report_A2 <- do.call(rbind, lapply(names(datasets_A2), function(plate_name) {

  A2 <- datasets_A2[[plate_name]]

  meta_sub <- meta_lookup[
    meta_lookup$Assay.Plate == plate_name &
    !is.na(meta_lookup$SampleName),
  ]

  # Extract quadrant from SampleName
  meta_sub$Quadrant <- regmatches(
    meta_sub$SampleName,
    regexpr("[ABCD]", meta_sub$SampleName)
  )

  quad_means <- sapply(c("A", "B", "C", "D"), function(q) {

    wells <- meta_sub$Array.Well[meta_sub$Quadrant == q]
    wells <- intersect(wells, colnames(A2))

    if (length(wells) == 0) return(NA)

    sample_counts <- colSums(A2[, wells, drop = FALSE] > 0, na.rm = TRUE)

    mean(sample_counts, na.rm = TRUE)
  })

  data.frame(
    Dataset    = plate_name,
    Quadrant_A = round(quad_means["A"], 2),
    Quadrant_B = round(quad_means["B"], 2),
    Quadrant_C = round(quad_means["C"], 2),
    Quadrant_D = round(quad_means["D"], 2),
    row.names  = NULL
  )
}))

print(report_A2)
write.csv(report_A2, "quadrant_avg_genes_report.csv", row.names = FALSE)

### Got Reads by multiplying reads/UMI by UMIs detected 
## Calculated Read/Gene Detected using that information 

## QC for all plates reported 
