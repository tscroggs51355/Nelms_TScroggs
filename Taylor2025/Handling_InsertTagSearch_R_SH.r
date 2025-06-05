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



### pUNOS Check - TS_March2025
filtered = "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data"

text_files_pUNos <- list.files(filtered, pattern = "_pUNos\\.txt$", full.names = TRUE)

data_list_pUNos <- lapply(text_files_pUNos, function(file) {
    sample_name <- gsub("_pUNos.txt", "", basename(file))

    # Check if file is empty
    if (file.info(file)$size > 0) {  
        df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence"))
        df$Sample <- sample_name
    } else {
        df <- data.frame(Count = NA, Sequence = NA, Sample = sample_name)  # Assign NA values
    }
    
    return(df)
})

final_data_pUNos <- do.call(rbind, data_list_pUNos)

# Save as CSV
write.csv(final_data_pUNos, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/pUNos_individual.csv", row.names = FALSE)

### pUbSplice Check - TS_March2025
filtered = "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data"

text_files_pUbSplice <- list.files(filtered, pattern = "_pUbSplice\\.txt$", full.names = TRUE)

data_list_pUbSplice <- lapply(text_files_pUbSplice, function(file) {
    sample_name <- gsub("_pUNos.txt", "", basename(file))

    # Check if file is empty
    if (file.info(file)$size > 0) {  
        df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence"))
        df$Sample <- sample_name
    } else {
        df <- data.frame(Count = NA, Sequence = NA, Sample = sample_name)  # Assign NA values
    }
    
    return(df)
})

final_data_pUbSplice <- do.call(rbind, data_list_pUbSplice)

# Save as CSV
write.csv(final_data_pUbSplice, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/pUbSplice_individual.csv", row.names = FALSE)

### NTS1 
text_files <- list.files("C:/Users/taylo/Desktop/June2025_Sequencing/NTS1_NTS1_1", pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(text_files, function(file) {
    sample_name <- sub("\\.txt$", "", basename(file))  

    if (file.info(file)$size > 0) {  
        df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence"))
        df$Sample <- sample_name
    } else {
        df <- data.frame(Count = NA, Sequence = NA, Sample = sample_name)  
    }

    return(df)  
})  
final_data <- do.call(rbind, data_list)  

write.csv(final_data, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG_NTS1_NTS1_1.csv", row.names = FALSE)


### SingleTubeNTS1-NTS1-1 

C:\Users\taylo\Desktop\June2025_Sequencing\singletube_NTS1-NTS1-1
text_files <- list.files("C:/Users/taylo/Desktop/June2025_Sequencing/singletube_NTS1-NTS1-1", pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(text_files, function(file) {
    sample_name <- sub("\\.txt$", "", basename(file))  

    if (file.info(file)$size > 0) {  
        df <- read.table(file, header = FALSE, col.names = c("Count", "Sequence"))
        df$Sample <- sample_name
    } else {
        df <- data.frame(Count = NA, Sequence = NA, Sample = sample_name)  
    }

    return(df)  
})  
final_data <- do.call(rbind, data_list)  

write.csv(final_data, "C:/Users/taylo/Desktop/TS_March2025/Mapped_Data/GGGTGGGCGCG_singtube_NTS1_NTS1_1.csv", row.names = FALSE)

