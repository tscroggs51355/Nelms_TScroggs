### TF Array 

setwd("C:/Users/taylo/Desktop/TF Project Master Doucments/Array") 
TFome = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Maize_TFome_Bulk_data _ corrected_BDN.csv")
GeneID = read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/V5GeneID.csv")
all_batches <- read.csv("all_batches_Locations.csv")

library(dplyr)
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]
library(stringr)
NTS1_Array <- read.csv("NTS1_Array.csv")
NTS2_Array <- read.csv("NTS2_Array.csv")
NTS3_Array <- read.csv("NTS3_Array.csv")
NTS4_Array <- read.csv("NTS4_Array.csv")
NTS5_Array <- read.csv("NTS5_Array.csv")
NTS6_Array <- read.csv("NTS6_Array.csv")
NTS7_Array <- read.csv("NTS7_Array.csv")
NTS8_Array <- read.csv("NTS8_Array.csv")

#sapply(all_locations, function(x) any(grepl(" ", x))) ## Looking for trailing white space 
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


plate12 = data.frame(array1 = c(rep(paste('NTS',1:5,sep=''),each = 2), rep(paste('NTS',6:13,sep=''),each = 2)), array2 = c(rep(paste('NTS',1:5,sep=''),each = 2)[c(2:10,1)], rep(paste('NTS',6:13,sep=''),each = 2)[c(2:16,1)]))
plate12$assayplate = apply(plate12, 1, paste, collapse = '-')

plate12 = plate12[(plate12$array1 %in% all_locations$Array.Plate) & (plate12$array2 %in% all_locations$Array.Plate),]

assay1 = lapply(plate12$array1, function(pl) { all_locations[all_locations[,2] == pl,c(3,1)]})
table(unlist(lapply(assay1,nrow))) #make sure it is all 176
assay1 = cbind(rep(plate12$assayplate, each = 176), do.call(rbind,assay1))
colnames(assay1)[1] = 'Assay.Plate'

assay2 = lapply(plate12$array2, function(pl) {
    out = all_locations[all_locations[,2] == pl,c(3,1)]
    out[,1] = LayoutPlate[out[,1]]
    return(out)
    })
table(unlist(lapply(assay2,nrow))) #make sure it is all 176
assay2 = cbind(rep(plate12$assayplate, each = 176), do.call(rbind,assay2))
colnames(assay2)[1] = 'Assay.Plate'

assay = rbind(assay1,assay2)
assay = assay[order(assay[,1],assay[,2]),]

      write.csv(assay, "AssayPlate_MetadataFull.csv", row.names = F)

metadata = assay 

metadata$TFomeStockID <- all_batches$TFomeStockID[
  match(metadata$SampleID, all_batches$SampleID)
]

V5ID <- GeneID$gene.ID[
  match(metadata$TFomeStockID, GeneID$clone)
]

TFomeID <- TFome$Gene.model[
  match(metadata$TFomeStockID, TFome$Stock.number)
]

metadata$Gene.model <- ifelse(
  is.na(V5ID),
  TFomeID,
  V5ID
)

write.csv(metadata, ""C:/Users/taylo/Desktop/TF Project Master Doucments/Array/metadata_TFAtlas.csv", row.names = FALSE)

## Metadata Updated for up to NTS8, but doesn't have any intermixed plate 

######################### Synbio Batch Information 
setwd("C:/Users/taylo/Desktop/TF Project Master Doucments") 

Batch1_Synbio <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Batch1_Synbio_FinalLocations.csv")
Batch2_Synbio <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Batch2_Synbio.csv")
Batch3_Synbio <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Batch3_Synbio.csv")
Batch4_Synbio <- read.csv("C:/Users/taylo/Desktop/TF Project Master Doucments/Batch4_Synbio_NTS5.csv")

all_batches <- bind_rows(
  Batch1_Synbio %>% select(SampleID, TFomeStockID),
  Batch2_Synbio %>% select(SampleID, TFomeStockID),
  Batch3_Synbio %>% select(SampleID, TFomeStockID),
  Batch4_Synbio %>% select(SampleID, TFomeStockID)
)

all_batches <- subset(all_batches, SampleID != "" & TFomeStockID != "")

sum(is.na(all_batches))

dim(all_batches)

all_batches_locations <- all_batches %>%
  distinct() %>%   
  left_join(
    all_locations %>% distinct(), 
    by = "SampleID",
  )


      write.csv(all_batches_Locations, "all_batches_Locations.csv", row.names = FALSE)


sum(all_batches_Locations %>%
  distinct() %>%                 
  count(SampleID))

### 2000 Total TFs 

5 Doubles in Batch 1 -- 5 samples in NTS1 and NTS2 Array Plates 

all_batches_Locations %>%
  count(SampleID) %>%
  filter(n > 1)

       SampleID n
1 S031591-137 4
2 S031591-139 4
3 S031591-147 4
4 S031591-162 4
5 S031591-168 4













setwd("C:/Users/taylo/Desktop/TF Project Master Doucments")
Batch2_Sequence <- read.csv("S036682_Plasmidpreparationmaster.csv")
colnames(Batch2_Sequence)[colnames(Batch2_Sequence) == "Sample.ID"] <- "SampleID"

all_batches <- left_join(
  all_batches,
  Batch2_Sequence %>% select(SampleID, Sequence),
  by = "SampleID"
)

FullCloningList <- left_join(
  all_batches,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)

dup_ids <- unique(all_batches$TFomeStockID[duplicated(all_batches$TFomeStockID)])
dup_ids <- unique(FullCloningList$SampleID[duplicated(FullCloningList$SampleID)])


dup_rows <- FullCloningList[FullCloningList$SampleID %in% dup_ids, ]

dup_rows ## Ask Brad 







      write.csv(FullCloningList, "FullCloningList.csv", row.names = FALSE)


























### NTS5 TF Array 
setwd("C:/Users/taylo/Desktop/2025_Nelms/August 2025") 
NTS5_Array <- read.csv("NTS5_Array_V2.csv")
Batch2_Synbio <- read.csv("Batch2_Synbio_NTS4.csv")
Batch3_Synbio <- read.csv("Batch3_Synbio_NTS4.csv")
Batch4_Synbio <- read.csv("Batch4_Synbio_NTS4.csv")

library(dplyr)
library(stringr)

all_locations <- bind_rows(
  NTS2_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS3_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS4_Array_Clean %>% select(SampleID, Library.Plate, Column, Row), 
  NTS5_Array %>% select(SampleID, Library.Plate, Column, Row))


Batch2_Synbio <- left_join(
  Batch2_Synbio,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)

      write.csv(Batch2_Synbio, "Batch2_Synbio_NTS5.csv", row.names = FALSE)


Batch3_Synbio <- left_join(
  Batch3_Synbio,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)
      write.csv(Batch3_Synbio, "Batch3_Synbio_NTS5.csv", row.names = FALSE)

Batch4_Synbio <-left_join(
  Batch4_Synbio,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)
      write.csv(Batch4_Synbio, "Batch4_Synbio_NTS5.csv", row.names = FALSE)
















###############################################################################################################################################################################################
## Directory for all files 
setwd("C:/Users/taylo/Desktop/2025_Nelms/August 2025") 

## NTS4 TF Information 
library(dplyr)
library(stringr)

NTS4_Array <- read.csv("NTS4_Array_clean.csv")
names(NTS4_Array)[names(NTS4_Array) == "Synbio.ID.Batch"] <- "SampleID"

names(Batch4_Synbio)[names(Batch4_Synbio) == "Project.ID"] <- "SampleID"
names(Batch4_Synbio)[names(Batch4_Synbio) == "Gene.name"] <- "TFomeStockID"


Batch2_Synbio <- read.csv("Batch2_Synbio_NTS3.csv")
Batch3_Synbio <- read.csv("Batch3_Synbio_NTS3.csv")
Batch4_Synbio <- read.csv("Batch4_Synbio_NTS4.csv")

names(Batch4_Synbio)[names(Batch4_Synbio) == "Project.ID"] <- "SampleID"
names(Batch4_Synbio)[names(Batch4_Synbio) == "TFomeStockID.x"] <- "TFomeStockID"

combined_synbio <- bind_rows(Batch2_Synbio, Batch3_Synbio, Batch4_Synbio) %>% 
select(SampleID, TFomeStockID)

combined_synbio_clean <- combined_synbio %>%
  filter(SampleID != "") %>%
  distinct(SampleID, .keep_all = TRUE)

NTS4_Array_Clean <- NTS4_Array %>%
  left_join(combined_synbio_clean, by = "SampleID")
## write.csv(NTS4_Array_Clean, "NTS4_Array_Clean.csv", row.names = FALSE)



colnames(NTS3_Array)[colnames(NTS3_Array) == "TFomeStockID.x"] <- "TFomeStockID"
colnames(Batch2_Synbio)[colnames(Batch2_Synbio) == "Library.plate"] <- "Library.Plate"
colnames(Batch3_Synbio)[colnames(Batch3_Synbio) == "Library.plate"] <- "Library.Plate"



setwd("C:/Users/taylo/Desktop/2025_Nelms/June 2025") 
NTS1_Array <- read.csv("NTS1_Array_Cleaned.csv")
NTS2_Array <- read.csv("NTS2_Array_clean.csv")
NTS3_Array <- read.csv("NTS3_Array_clean.csv")



all_locations <- bind_rows(
  NTS2_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS3_Array %>% select(SampleID, Library.Plate, Column, Row), 
  NTS4_Array_Clean %>% select(SampleID, Library.Plate, Column, Row))



Batch1_Synbio <- Batch1_Synbio %>%
  left_join(all_locations, by = "SampleID")
  write.csv(Batch1_Synbio, "Batch1_Synbio_NTS2.csv", row.names = FALSE)

Batch2_Synbio <- left_join(
  Batch2_Synbio,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)

      write.csv(Batch2_Synbio, "Batch2_Synbio_NTS4.csv", row.names = FALSE)


Batch3_Synbio <- left_join(
  Batch3_Synbio,
  all_locations,
  by = "SampleID",
  suffix = c("", "")
)
      write.csv(Batch3_Synbio, "Batch3_Synbio_NTS4.csv", row.names = FALSE)

Batch4_Synbio <- Batch4_Synbio[, -((ncol(Batch4_Synbio) - 5):ncol(Batch4_Synbio))]
Batch4_Synbio <- Batch4_Synbio %>%
  left_join(NTS4_Array_Clean, by = "SampleID")
      write.csv(Batch4_Synbio, "Batch4_Synbio_NTS4.csv", row.names = FALSE)




















############################################################################################################################

## NTS2 TF Information 
NTS2_Array <- read.csv("NTS2_Array.csv")
NTS2_Array_clean <- NTS2_Array[, 1:5]
names(NTS2_Array_clean)[names(NTS2_Array_clean) == "Synbio.ID.Batch"] <- "SampleID"

NTS2_Array_clean <- NTS2_Array_clean[,c(4,1,2,3)]

head(NTS2_Array_clean)
tail(NTS2_Array_clean)


## NTS1 TF Information 
NTS1_Array <- read.csv("September2024_NTS1.csv")
NTS1_Array_clean <- NTS1_Array[,c(1,2,14,15)]
NTS1_Array_clean <- NTS1_Array_clean[!is.na(NTS1_Array_clean$Library.plate.location), ]
names(NTS1_Array_clean)[names(NTS1_Array_clean) == "ProjectID"] <- "SampleID"
names(NTS1_Array_clean)[names(NTS1_Array_clean) == "Stock.number"] <- "TFomeStockID"
NTS1_Array_clean$Row <- sub("([A-Z]+)([0-9]+)", "\\1", NTS1_Array_clean$Library.plate.location)
NTS1_Array_clean$Column <- sub("([A-Z]+)([0-9]+)", "\\2", NTS1_Array_clean$Library.plate.location)
NTS1_Array_clean$Column <- as.numeric(NTS1_Array_clean$Column)
NTS1_Array_clean <- NTS1_Array_clean[,c(1,2,3,6,5)]

library(dplyr)
# mCherry Samples Added 
mCherry_rows <- NTS2_Array_clean %>%
  filter(SampleID == "mCherry") %>%
  rename(Library.plate = Library.Plate) %>%
  mutate(Library.plate = "NTS1")
  NTS1_Array_clean <- bind_rows(NTS1_Array_clean, mCherry_rows)
head(NTS1_Array_clean)
tail(NTS1_Array_clean)

write.csv(NTS1_Array_clean, "NTS1_Array_Cleaned.csv", row.names = FALSE)



## Batch 1 Synbio
Batch1_Synbio <- read.csv("S031591list.csv")

library(dplyr)
library(stringr)

Batch1_Synbio <- Batch1_Synbio %>%
  mutate(TFomeStockID = str_extract(COA.File.Name, "(?<=-)[^-.]+(?=\\.zip)"))

  Batch1_Synbio <- Batch1_Synbio[,c(1,3)]
names(Batch1_Synbio)[names(Batch1_Synbio) == "ProjectID"] <- "SampleID"

## Batch 2 Synbio 
Batch2_Synbio <- read.csv("S036682_Batch2_master.csv")
Batch2_Synbio <- Batch2_Synbio[,1:3]
names(Batch2_Synbio)[names(Batch2_Synbio) == "Sample.ID"] <- "SampleID"
names(Batch2_Synbio)[names(Batch2_Synbio) == "gene.Name"] <- "TFomeStockID"

#Batch2_Synbio_Sub2 <- read.csv("S036682_Batch2_Sub2.csv")
#Batch2_Synbio_Sub2 <- Batch2_Synbio_Sub2[,1:4]
#names(Batch2_Synbio_Sub2)[names(Batch2_Synbio_Sub2) == "Sample.ID"] <- "SampleID"
#names(Batch2_Synbio_Sub2)[names(Batch2_Synbio_Sub2) == "gene.Name"] <- "TFomeStockID"


## Batch 3 Synbio 
Batch3_Synbio <- read.csv("S034674_Batch3.csv")
Batch3_Synbio <- Batch3_Synbio[,1:4]
names(Batch3_Synbio)[names(Batch3_Synbio) == "Sample.ID"] <- "SampleID"
names(Batch3_Synbio)[names(Batch3_Synbio) == "Gene.Name"] <- "TFomeStockID"

library(dplyr)
combined_synbio <- bind_rows(Batch1_Synbio, Batch2_Synbio, Batch2_Synbio_Sub2,Batch3_Synbio) %>%
  select(SampleID, TFomeStockID)

combined_synbio_clean <- combined_synbio %>%
  filter(SampleID != "") %>%
  distinct(SampleID, .keep_all = TRUE)

NTS2_Array_clean <- NTS2_Array_clean %>%
  left_join(combined_synbio_clean, by = "SampleID")
write.csv(NTS2_Array_clean, "NTS2_Array_clean.csv", row.names = FALSE)


all_locations <- bind_rows(
  NTS1_Array_clean %>% select(SampleID, Library.plate, Column, Row),
  NTS2_Array_clean %>% rename(Library.plate = Library.Plate) %>%
    select(SampleID, Library.plate, Column, Row)
)

Batch1_Synbio <- Batch1_Synbio %>%
  left_join(all_locations, by = "SampleID")
  write.csv(Batch1_Synbio, "Batch1_Synbio_NTS2.csv", row.names = FALSE)

Batch2_Synbio <- Batch2_Synbio %>%
  left_join(all_locations, by = "SampleID")
      write.csv(Batch2_Synbio, "Batch2_Synbio_NTS2.csv", row.names = FALSE)


Batch3_Synbio <- Batch3_Synbio %>%
  left_join(all_locations, by = "SampleID")
      write.csv(Batch3_Synbio, "Batch3_Synbio_NTS2.csv", row.names = FALSE)
