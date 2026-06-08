## Final Array 


setwd("C:/Users/taylo/Desktop/TF Project Master Doucments") 

metadata = read.csv("Metadata_FullPlate_Jan2026.csv")

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/Array") 
NTS1_Array <- read.csv("NTS1_Array.csv")
NTS2_Array <- read.csv("NTS2_Array.csv")
NTS3_Array <- read.csv("NTS3_Array.csv")
NTS4_Array <- read.csv("NTS4_Array.csv")
NTS5_Array <- read.csv("NTS5_Array.csv")
NTS6_Array <- read.csv("NTS6_Array.csv")
NTS7_Array <- read.csv("NTS7_Array.csv")
NTS8_Array <- read.csv("NTS8_Array.csv")
library(dplyr)

colnames(NTS1_Array)[colnames(NTS1_Array) == "Library.plate"] <- "Library.Plate"
NTS1_Array <- NTS1_Array[order(NTS1_Array$Column, NTS1_Array$Row), ]

all_locations <- bind_rows(
  NTS1_Array %>% select(SampleID, Library.Plate, Column, Row),
  NTS2_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS3_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS4_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS5_Array %>% select(SampleID, Library.Plate, Column, Row),
  NTS6_Array %>% select(SampleID, Library.Plate, Column, Row),
  NTS7_Array %>% select(SampleID, Library.Plate, Column, Row),
  NTS8_Array %>% select(SampleID, Library.Plate, Column, Row),
  )

all_locations <- dplyr::rename(all_locations, Array.Plate = Library.Plate)
all_locations <- dplyr::rename(all_locations, Array.Column = Column)
all_locations <- dplyr::rename(all_locations, Array.Row = Row)

all_locations$Array.Well <- paste0(
  all_locations$Array.Row,
  sprintf("%02d", all_locations$Array.Column)
)

all_locations[] <- lapply(all_locations, function(x) {
    if (is.character(x)) trimws(x) else x
})

all_locations <- all_locations[, !colnames(all_locations) %in% c("Array.Column", "Array.Row")]



SynbioFiles = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/SynbioFiles/FullPlasmidListv0.5.csv")

## --- Build location string in all_locations ---
all_locations$Location <- paste(all_locations$Array.Plate,
                                all_locations$Array.Well,
                                sep = "-")

## --- For each plasmid, collect all matching locations ---
pl <- SynbioFiles$Plasmid

loc_list <- lapply(pl, function(x) {
  idx <- which(all_locations$SampleID == x)
  if (length(idx) == 0) return(NA_character_)
  all_locations$Location[idx]
})

## --- Determine max number of matches ---
max_locs <- max(sapply(loc_list, length))

## --- Expand list into a matrix with padded NAs ---
loc_mat <- do.call(rbind, lapply(loc_list, function(v) {
  length(v) <- max_locs
  v
}))

colnames(loc_mat) <- paste0("Location_", seq_len(max_locs))

## --- Bind to SynbioFiles 
SynbioFiles_with_locations <- cbind(SynbioFiles, loc_mat)

#### Counting how many rows where both location_1 and Location_2 are empty, this would indicate they have not been arrayed 
> sum(is.na(SynbioFiles_with_locations$Location_1) &
+     is.na(SynbioFiles_with_locations$Location_2))
[1] 804
> 
> table(all_locations$Array.Plate) ## How many of each plate are there 

NTS1 NTS2 NTS3 NTS4 NTS5 NTS6 NTS7 NTS8 
 176  176  176  176  176  176  176  176 


## > unique(all_locations$SampleID[duplicated(all_locations$SampleID)])
[1] "mCherry"     "S031591-137" "S031591-139" "S031591-147" "S031591-162"
[6] "S031591-168"
> 


## Subsetting just to NTS1-NTS5 
all_locations_sub <- all_locations[
  all_locations$Array.Plate %in% paste0("NTS", 1:5),
]


all_locations_sub <- all_locations_sub[ all_locations_sub$SampleID != "mCherry", ]

unique(all_locations_sub$SampleID[duplicated(all_locations_sub$SampleID)])

rep_ids = unique(all_locations_sub$SampleID[duplicated(all_locations_sub$SampleID)])




## all SampleIDs in NTS1–NTS5
sub_ids <- all_locations_sub$SampleID

## subtract replicated ones
unique_nonrep_ids <- setdiff(sub_ids, rep_ids)

## optional: show them
unique_nonrep_ids ## List od IDs to be randomly selected from for the remaining spots  


###### 

#write.csv(SynbioFiles_with_locations, "SynbioFiles_with_locations.csv")
#write.csv(all_locations, "ArrayLocations_NTS1_NTS8.csv")

all_locations = NTS1, NTS8 Array (1408 Total including Mcherry) (1307 if substracting the mCherry and repeats) 
all_locations_sub = NTS1, NTS5 Array (820 minus the mCherry) 
SynbioFiles_with_locations = AllSynbio (5 things replicated with two locations on the array plate)

SynbioFiles_with_locations
None = 804 
Location 1 = 1303 
Location 1 + 2 = 5 
Location 2 = 0 


Total TFs = 2112 

## One sample with duplicates in SynbioFiles_with_locations, S036682-515 appears twice, removed these in Synbio_no_dupes
none  <- sum(is.na(Location_1 ocation_2))
only1 <- sum(!is.na(Location_1) &  is.na(Location_2))
only2 <- sum( is.na(Location_1) & !is.na(Location_2))
both  <- sum(!is.na(Location_1) & !is.na(Location_2))

none; only1; only2; both
none + only1 + only2 + both

Total for 13 plates: 2132 
Currently Arrayed: 1307 
To be Arrayed: 804 

## Randomly Sampling the Subset of NTS1, NTS5 Above 
set.seed(1)   # optional, only if you want reproducibility
subsample_30 <- sample(unique_nonrep_ids, 30)

Synbio_no_dupes <- SynbioFiles_with_locations[
    !duplicated(SynbioFiles_with_locations$Plasmid),
]

ids_none <- Synbio_no_dupes$Plasmid[
    is.na(Synbio_no_dupes$Location_1) &
    is.na(Synbio_no_dupes$Location_2)
]

RemainingTFS = sample(c(ids_none,sample(subsample_30,16)))

BlankArray = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Blank_Array.csv")

### 1. Identify empty wells in BlankArray
empty_wells <- which(is.na(BlankArray$SampleID) | BlankArray$SampleID == "")

n_per_plate <- length(empty_wells)

if (length(RemainingTFS) < 5 * n_per_plate) {
    stop("Not enough TFs to fill 5 plates.")
}

set.seed(123)
shuffled_TFs <- sample(RemainingTFS, size = 5 * n_per_plate, replace = FALSE)

TF_chunks <- split(shuffled_TFs, rep(1:5, each = n_per_plate))

fill_plate <- function(blank, tf_list, plate_name) {
    df <- blank
    df$Library.Plate <- plate_name
    df$SampleID[empty_wells] <- tf_list
    df$TFomeStockID[empty_wells] <- NA
    df
}

NTS9  <- fill_plate(BlankArray, TF_chunks[[1]], "NTS9")
NTS10 <- fill_plate(BlankArray, TF_chunks[[2]], "NTS10")
NTS11 <- fill_plate(BlankArray, TF_chunks[[3]], "NTS11")
NTS12 <- fill_plate(BlankArray, TF_chunks[[4]], "NTS12")
NTS13 <- fill_plate(BlankArray, TF_chunks[[5]], "NTS13")

all_TFs <- c(
    NTS9$SampleID, NTS10$SampleID, NTS11$SampleID,
    NTS12$SampleID, NTS13$SampleID
)

stopifnot(
    length(unique(all_TFs[all_TFs != "mCherry"])) ==
    length(all_TFs[all_TFs != "mCherry"])
)

### ------------------------------------------------------------
### 0.  Build lookup tables from the five new plates
### ------------------------------------------------------------

make_lookup <- function(df) {
    df$Well <- paste0(df$Row, sprintf("%02d", df$Column))
    df$Location <- paste0(df$Library.Plate, "-", df$Well)
    df[, c("SampleID", "Location")]
}

lookup_NTS9  <- make_lookup(NTS9)
lookup_NTS10 <- make_lookup(NTS10)
lookup_NTS11 <- make_lookup(NTS11)
lookup_NTS12 <- make_lookup(NTS12)
lookup_NTS13 <- make_lookup(NTS13)

lookup_all <- rbind(
    lookup_NTS9,
    lookup_NTS10,
    lookup_NTS11,
    lookup_NTS12,
    lookup_NTS13
)

lookup_all <- lookup_all[lookup_all$SampleID != "mCherry", ]


### ------------------------------------------------------------
### 1. Merge lookup into SynbioFiles_with_locations
### ------------------------------------------------------------

SynbioFiles_with_locations <- merge(
    SynbioFiles_with_locations,
    lookup_all,
    by.x = "Plasmid",
    by.y = "SampleID",
    all.x = TRUE
)

### ------------------------------------------------------------
### 2. Fill Location_1 only where missing
### ------------------------------------------------------------

### Case 1: Location_1 is empty → fill Location_1
fill_L1 <- is.na(SynbioFiles_with_locations$Location_1) &
           !is.na(SynbioFiles_with_locations$Location)

SynbioFiles_with_locations$Location_1[fill_L1] <-
    SynbioFiles_with_locations$Location[fill_L1]

write.csv(NTS9,  "NTS9_Array_May26.csv",  row.names = FALSE)
write.csv(NTS10, "NTS10_Array_May26.csv", row.names = FALSE)
write.csv(NTS11, "NTS11_Array_May26.csv", row.names = FALSE)
write.csv(NTS12, "NTS12_Array_May26.csv", row.names = FALSE)
write.csv(NTS13, "NTS13_Array_May26.csv", row.names = FALSE)

write.csv(SynbioFiles_with_locations, "SynbioSequencing_Locations_filled_May2026.csv", row.names = FALSE)

write.csv(subsample_30, "SubsampledNTS1,NTS5.csv", row.names = FALSE)