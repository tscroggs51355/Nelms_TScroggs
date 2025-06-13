setwd("C:/Users/taylo/Desktop")

Batch2_Synbio <- read.csv("Batch2_Synbio.csv")
Batch2_Synbio <- Batch2_Synbio %>%
  separate(Sample.ID, into = c("Batch", "SampleNumber"), sep = "-")
  
Batch2_Synbio <- na.omit(Batch2_Synbio)
head(Batch2_Synbio)
write.csv(Batch2_Synbio, "Batch2_Synbio.csv")


Batch3_Synbio <- read.csv("Batch3_Synbio.csv")
library(tidyr)
Batch3_Synbio <- Batch3_Synbio %>%
  separate(Sample.ID, into = c("Batch", "SampleNumber"), sep = "-")
  
Batch3_Synbio <- na.omit(Batch3_Synbio)
head(Batch3_Synbio)
write.csv(Batch3_Synbio, "Batch3_Synbio.csv")


