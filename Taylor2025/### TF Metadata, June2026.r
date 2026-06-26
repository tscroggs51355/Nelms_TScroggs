### TF Metadata, June2026

TFome = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Maize_TFome_Bulk_data _ corrected_BDN.csv")
GeneID = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/V5GeneID.csv")
metadata = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Metadata_FullPlate_Jan2026.csv")

plate_rows <- LETTERS[1:16]          # A:P
plate_cols <- sprintf("%02d", 2:23)  # 02:23

wells <- as.vector(outer(plate_rows, plate_cols, paste0))

# Quadrant:
# A:B = A, C:D = B, E:F = C, G:H = D,
# I:J = A, K:L = B, M:N = C, O:P = D
quadrant_by_row <- rep(c("A", "B", "C", "D"), times = 2, each = 2)

# Barcode position:
# rows A:B = barcode group 1
# rows C:D = barcode group 2
# ...
# rows O:P = barcode group 8
barcode_group_by_row <- rep(1:8, each = 2)

# Columns 02:12 and 13:23 are two matching barcode blocks
barcode_number_by_col <- rep(0:10, times = 2)

# Build lookup for every well
meta <- expand.grid(
  row = plate_rows,
  col = plate_cols,
  stringsAsFactors = FALSE
)

meta$Array.Well <- paste0(meta$row, meta$col)

meta$Quadrant <- quadrant_by_row[match(meta$row, plate_rows)]

meta$BC <- paste0(
  9 + barcode_group_by_row[match(meta$row, plate_rows)] - 1 +
    barcode_number_by_col[match(meta$col, plate_cols)] * 8,
  "s"
)

# Keep only needed columns
meta <- meta[, c("Array.Well", "BC", "Quadrant")]

# Add well information to your metadata
metadata$BC <- meta$BC[match(metadata$Array.Well, meta$Array.Well)]
metadata$Quadrant <- meta$Quadrant[match(metadata$Array.Well, meta$Array.Well)]

# Fill only missing SampleName values
missing_name <- is.na(metadata$SampleName) | metadata$SampleName == ""

metadata$SampleName[missing_name] <- paste(
  metadata$Assay.Plate[missing_name],
  metadata$Quadrant[missing_name],
  metadata$BC[missing_name],
  sep = "_"
)

# Optional: remove helper columns
metadata$BC <- NULL
metadata$Quadrant <- NULL

meta[meta$Array.Well %in% c("A02", "B02", "A13", "B13"), ]

write.csv(metadata, "metadata_updated_June2026.csv")
##Updated with information to the NTS8-NTS8 plates 
