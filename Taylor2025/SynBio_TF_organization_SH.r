## Directory for all files 
setwd("C:/Users/taylo/Desktop/2025_Nelms/June 2025")

## NTS2 TF Information 
NTS2_Array <- read.csv("NTS2_Array.csv")
NTS2_Array_clean <- NTS2_Array[, 1:5]
names(NTS2_Array_clean)[names(NTS2_Array_clean) == "Synbio.ID.Batch"] <- "SampleID"

NTS2_Array_clean <- NTS2_Array_clean[,c(4,1,2,3)]

head(NTS2_Array_clean)
tail(NTS2_Array_clean)

## NTS3 TF Information 
NTS3_Array_clean <- read.csv("NTS3_Array.csv")
names(NTS3_Array_clean)[names(NTS3_Array_clean) == "Synbio.ID.Batch"] <- "SampleID"

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
combined_synbio <- bind_rows(Batch1_Synbio, Batch2_Synbio,Batch3_Synbio) %>%
  select(SampleID, TFomeStockID)

combined_synbio_clean <- combined_synbio %>%
  filter(SampleID != "") %>%
  distinct(SampleID, .keep_all = TRUE)

NTS2_Array_clean <- NTS2_Array_clean %>%
  left_join(combined_synbio_clean, by = "SampleID")
write.csv(NTS2_Array_clean, "NTS2_Array_clean.csv", row.names = FALSE)

NTS3_Array_clean <- NTS3_Array_clean %>%
  left_join(combined_synbio_clean, by = "SampleID")
write.csv(NTS3_Array_clean, "NTS3_Array_clean.csv")

all_locations <- bind_rows(
  NTS1_Array_clean %>%
    select(SampleID, Library.plate, Column, Row),
  NTS2_Array_clean %>%
    rename(Library.plate = Library.Plate) %>%
    select(SampleID, Library.plate, Column, Row),
  NTS3_Array_clean %>%
    rename(Library.plate = Library.Plate) %>%
    select(SampleID, Library.plate, Column, Row)
)

Batch1_Synbio <- Batch1_Synbio %>%
  left_join(all_locations, by = "SampleID")
  write.csv(Batch1_Synbio, "Batch1_Synbio_NTS3.csv", row.names = FALSE)

Batch2_Synbio <- Batch2_Synbio %>%
  left_join(all_locations, by = "SampleID")
      write.csv(Batch2_Synbio, "Batch2_Synbio_NTS3.csv", row.names = FALSE)


Batch3_Synbio <- Batch3_Synbio %>%
  left_join(all_locations, by = "SampleID")
      write.csv(Batch3_Synbio, "Batch3_Synbio_NTS3.csv", row.names = FALSE)
