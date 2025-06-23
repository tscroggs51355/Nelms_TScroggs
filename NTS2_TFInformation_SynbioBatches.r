
setwd("C:/Users/taylo/Desktop/2025_Nelms/June 2025")

Batch2_Synbio <- read.csv("Batch2_Synbio.csv")
Batch3_Synbio <- read.csv("Batch3_Synbio.csv")
NTS2_Array <- read.csv("NTS2_Array.csv")


library(dplyr)

Batch2_Synbio_NTS2 <- Batch2_Synbio %>%
  left_join(NTS2_Array, by = c("Sample.ID" = "Synbio.ID.Batch")) %>%
  select(Sample.ID, everything()) 
Batch2_Synbio_NTS2 <- Batch2_Synbio_NTS2[,1:5]
head(Batch2_Synbio_NTS2)
write.csv(Batch2_Synbio_NTS2, file = "Batch2_Synbio_NTS2.csv", row.names = FALSE)



Batch3_Synbio_NTS2 <- Batch3_Synbio %>%
  left_join(NTS2_Array, by = c("Sample.ID" = "Synbio.ID.Batch")) %>%
  select(Sample.ID, everything()) 
Batch3_Synbio_NTS2 <- Batch3_Synbio_NTS2[,1:5]
head(Batch3_Synbio_NTS2)
write.csv(Batch3_Synbio_NTS2, file = "Batch3_Synbio_NTS2.csv", row.names = FALSE)
