library(tidyverse)

filtered = "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data"

text_files <- list.files("C:/Users/taylo/Desktop/TS_March2025/Mapped_Data", pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(text_files, function(file) {
    sample_name <- gsub(".txt", "", basename(file))  
    df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence")) 
    df$Sample <- sample_name 
    return(df)
})


final_data <- do.call(rbind, data_list)  

write.csv(final_data, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG_individual.csv", row.names = FALSE)


df_wide <- reshape(final_data, 
                   idvar = "Sequence", 
                   timevar = "Sample", 
                   direction = "wide")

colnames(df_wide) <- gsub("Count.", "", colnames(df_wide)) 

head(df_wide)

write.csv(df_wide, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG.csv", row.names = FALSE)
