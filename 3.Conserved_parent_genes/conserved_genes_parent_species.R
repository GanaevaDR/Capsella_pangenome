# Establish conserved genes from parent species


# process mosdepth output files. Merge files per sample into one table

# load libraries
library(readr)
library(dplyr)

rm(list = ls())
folder_path <- "C:/Users/Dasha/Desktop/masters/parent_species/gene_coverage_mosdepth/Cori"
file_list <- list.files(path = folder_path, pattern = "*.bed.gz", full.names = TRUE)

# Initialize an empty dataframe for merging
merged_df <- data.frame()

# Loop through each file to read and process
for (file in file_list) {
  df <- read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  # Extract the relevant columns
  temp_df <- df %>%
    select(4, 5) %>%  
    rename(key = X4, value = X5) 
  
  # If merged_df is empty, initialize it with the first dataframe
  if (nrow(merged_df) == 0) {
    merged_df <- temp_df
  } else {
    # Merge with existing dataframe by 'key'
    merged_df <- full_join(merged_df, temp_df, by = "key", suffix = c("", paste0("_", basename(file))))
  }
}

# Rename the value columns
colnames(merged_df)[-1] <- gsub("_mean.regions.bed.gz", "", list.files(path = folder_path, pattern = "*_mean.regions.bed.gz"))
colnames(merged_df)[1] <- "gene"

merged_vec <- unlist(merged_df[, -1])

Cori_df <- as.data.frame(merged_df)
rownames(Cori_df) <- Cori_df$gene
Cori_df$gene <- NULL

write.table(Cori_df, 'C_ori_mean_coverage.txt',
            sep = '\t', col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = 'NA')

#==============================================================
# C rubella
rm(list = ls())
# Median
folder_path <- "C:/Users/Dasha/Desktop/masters/parent_species/gene_coverage_mosdepth/Crub"
file_list <- list.files(path = folder_path, pattern = "*.bed.gz", full.names = TRUE)

# Initialize an empty dataframe for merging
merged_df <- data.frame()

# Loop through each file to read and process
for (file in file_list) {
  df <- read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  # Extract the relevant columns
  temp_df <- df %>%
    select(4, 5) %>% 
    rename(key = X4, value = X5)  
  
  if (nrow(merged_df) == 0) {
    merged_df <- temp_df
  } else {
    # Merge with existing dataframe by 'key'
    merged_df <- full_join(merged_df, temp_df, by = "key", suffix = c("", paste0("_", basename(file))))
  }
}

# Rename the value columns
colnames(merged_df)[-1] <- gsub("_mean.regions.bed.gz", "", list.files(path = folder_path, pattern = "*_mean.regions.bed.gz"))
colnames(merged_df)[1] <- "gene"

merged_vec <- unlist(merged_df[, -1])

Crub_df <- as.data.frame(merged_df)
rownames(Crub_df) <- Crub_df$gene
Crub_df$gene <- NULL

write.table(Crub_df, 'C_rub_mean_coverage.txt',
            sep = '\t', col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = 'NA')

#==========================================================================
# C. grandiflora

rm(list = ls())

folder_path <- "C:/Users/Dasha/Desktop/masters/parent_species/gene_coverage_mosdepth/Cgra"
file_list <- list.files(path = folder_path, pattern = "*.bed.gz", full.names = TRUE)


# Initialize an empty dataframe for merging
merged_df <- data.frame()

# Loop through each file to read and process
for (file in file_list) {
  # Read the current file
  df <- read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  # Extract the relevant columns
  temp_df <- df %>%
    select(4, 5) %>%  
    rename(key = X4, value = X5) 
  
  if (nrow(merged_df) == 0) {
    merged_df <- temp_df
  } else {
    # Merge with existing dataframe by 'key'
    merged_df <- full_join(merged_df, temp_df, by = "key", suffix = c("", paste0("_", basename(file))))
  }
}

# Rename the value columns
colnames(merged_df)[-1] <- gsub("_mean.regions.bed.gz", "", list.files(path = folder_path, pattern = "*_mean.regions.bed.gz"))
colnames(merged_df)[1] <- "gene"

merged_vec <- unlist(merged_df[, -1])

Cgra_df <- as.data.frame(merged_df)
rownames(Cgra_df) <- Cgra_df$gene
Cgra_df$gene <- NULL

write.table(Cgra_df, 'C_gra_mean_coverage.txt',
            sep = '\t', col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = 'NA')

#==========================================================================
#==========================================================================

# work with obtained tables

rm(list = ls())
library(dplyr)

Co_mean <- read.table("C_ori_mean_coverage.txt", header = TRUE, row.names = 1,
                      sep="\t", na.strings = "NA", dec = ".",
                      strip.white = TRUE)

Cr_mean <- read.table("C_rub_mean_coverage.txt", header = TRUE, row.names = 1,
                      sep="\t", na.strings = "NA", dec = ".",
                      strip.white = TRUE)

Cg_mean <- read.table("C_gra_mean_coverage.txt", header = TRUE, row.names = 1,
                      sep="\t", na.strings = "NA", dec = ".",
                      strip.white = TRUE)

# filter samples by gene coverage (depth)
Co_mean <- Co_mean %>%
  select(where(~ mean(. , na.rm = TRUE) >= 10)) #3 samples not passed

Cr_mean <- Cr_mean %>%
  select(where(~ mean(. , na.rm = TRUE) >= 10)) #4 samples not passed

Cg_mean <- Cg_mean %>%
  select(where(~ mean(. , na.rm = TRUE) >= 10)) #4 samples not passed

# calculate mean and threshold
Co_means <- apply(Co_mean, 2, mean)
Co_thr <- Co_means * 0.1

Cr_means <- apply(Cr_mean, 2, mean)
Cr_thr <- Cr_means * 0.1

Cg_means <- apply(Cg_mean, 2, mean)
Cg_thr <- Cg_means * 0.1

mean(Co_thr)
mean(Cr_thr)
mean(Cg_thr)

# filter out all values less than 10% of mean
Co_mean_clone <- Co_mean
Co_mean_clone$count <- apply(Co_mean, 1, function(row) sum(row > Co_thr))

Cr_mean_clone <- Cr_mean
Cr_mean_clone$count <- apply(Cr_mean, 1, function(row) sum(row > Cr_thr))

Cg_mean_clone <- Cg_mean
Cg_mean_clone$count <- apply(Cg_mean, 1, function(row) sum(row > Cg_thr))

# merge Cg and Cr 
Cr_Cg <- cbind(Cr_mean, Cg_mean)

Cr_Cg_means <- apply(Cr_Cg, 2, mean)
Cr_Cg_thr <- Cr_Cg_means * 0.1

Cr_Cg_mean_clone <- Cr_Cg
Cr_Cg_mean_clone$count <- apply(Cr_Cg, 1, function(row) sum(row > Cr_Cg_thr))

library(dplyr)
Co_filt <- filter(Co_mean_clone, count==ncol(Co_mean))
Co_filt$Co <- rownames(Co_filt)

Cr_Cg_filt <- filter(Cr_Cg_mean_clone, count==ncol(Cr_Cg))
Cr_Cg_filt$Cr <- rownames(Cr_Cg_filt)


library(readr)
orthopairs <- read_delim("orthopairs.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

orthopairs$Co <-  substr(orthopairs$Co, 1, nchar(orthopairs$Co) - 3)
orthopairs$Cr <-  substr(orthopairs$Cr, 1, nchar(orthopairs$Cr) - 3)

merged <- merge(orthopairs, Co_filt, by="Co")
merged1 <- merge(merged, Cr_Cg_filt, by="Cr")

merged <- merge(orthopairs, Cr_Cg_filt, by="Cr")
merged1 <- merge(merged, Co_filt, by="Co") #20,650 genes
merged1$count.x <- NULL
merged1$count.y <- NULL

write.table(merged1, 'conserved_genes_filt.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')
