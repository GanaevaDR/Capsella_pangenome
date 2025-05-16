library(readr)
library(dplyr)
library(stringr)

# Eng population

dam <- read_delim("Eng_damaged_exons.txt", 
                  delim = "\t", escape_double = FALSE, col_names = FALSE,
                  trim_ws = TRUE)

dam$Ath <- str_extract(dam$X1, "AT[^_]*")
dam$X11 <- gsub("\\.", " ", dam$Ath)
dam_vec <- unlist(strsplit(dam$X11, " "))
dam_d <- as.data.frame(dam_vec)
dam_d  <- na.omit(dam_d)

write.table(dam_d, 'fg_Eng_damaging.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE, na = 'NA')


GO <- read_delim("GO_Eng_damaging.txt", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)

colnames(GO) <- gsub(" ", "_", colnames(GO))

GO_f <- filter(GO, Fold_Enrichment >=2 & FDR< 0.05)

GO_f <- GO_f[order(GO_f$Fold_Enrichment, decreasing = TRUE),]  

write.table(GO_f, 'res_Eng_damaging.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')

#===============================================================================

# KBG population

dam <- read_delim("KBG_damaged_exons.txt", 
                  delim = "\t", escape_double = FALSE, col_names = FALSE,
                  trim_ws = TRUE)

dam$Ath <- str_extract(dam$X1, "AT[^_]*")
dam$X11 <- gsub("\\.", " ", dam$Ath)
dam_vec <- unlist(strsplit(dam$X11, " "))
dam_d <- as.data.frame(dam_vec)
dam_d  <- na.omit(dam_d)

write.table(dam_d, 'fg_KBG_damaging.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE, na = 'NA')


GO <- read_delim("GO_KBG_damaging.txt", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)

colnames(GO) <- gsub(" ", "_", colnames(GO))

GO_f <- filter(GO, Fold_Enrichment >=2 & FDR< 0.05)

GO_f <- GO_f[order(GO_f$Fold_Enrichment, decreasing = TRUE),]  

write.table(GO_f, 'res_KBG_damaging.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')


#===============================================================================

# le3 population

dam <- read_delim("le3_damaged_exons.txt", 
                  delim = "\t", escape_double = FALSE, col_names = FALSE,
                  trim_ws = TRUE)

dam$Ath <- str_extract(dam$X1, "AT[^_]*")
dam$X11 <- gsub("\\.", " ", dam$Ath)
dam_vec <- unlist(strsplit(dam$X11, " "))
dam_d <- as.data.frame(dam_vec)
dam_d  <- na.omit(dam_d)

write.table(dam_d, 'fg_le3_damaging.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE, na = 'NA')


GO <- read_delim("GO_le3_damaging.txt", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)

colnames(GO) <- gsub(" ", "_", colnames(GO))

GO_f <- filter(GO, Fold_Enrichment >=2 & FDR< 0.05)

GO_f <- GO_f[order(GO_f$Fold_Enrichment, decreasing = TRUE),]  

write.table(GO_f, 'res_le3_damaging.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')

#===============================================================================

# MSK population

dam <- read_delim("MSK_damaged_exons.txt", 
                  delim = "\t", escape_double = FALSE, col_names = FALSE,
                  trim_ws = TRUE)

dam$Ath <- str_extract(dam$X1, "AT[^_]*")
dam$X11 <- gsub("\\.", " ", dam$Ath)
dam_vec <- unlist(strsplit(dam$X11, " "))
dam_d <- as.data.frame(dam_vec)
dam_d  <- na.omit(dam_d)

write.table(dam_d, 'fg_MSK_damaging.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE, na = 'NA')


GO <- read_delim("GO_MSK_damaging.txt", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)

colnames(GO) <- gsub(" ", "_", colnames(GO))

GO_f <- filter(GO, Fold_Enrichment >=2 & FDR< 0.05)

GO_f <- GO_f[order(GO_f$Fold_Enrichment, decreasing = TRUE),]  

write.table(GO_f, 'res_MSK_damaging.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')

#===============================================================================

# Mur population

dam <- read_delim("Mur_damaged_exons.txt", 
                  delim = "\t", escape_double = FALSE, col_names = FALSE,
                  trim_ws = TRUE)

dam$Ath <- str_extract(dam$X1, "AT[^_]*")
dam$X11 <- gsub("\\.", " ", dam$Ath)
dam_vec <- unlist(strsplit(dam$X11, " "))
dam_d <- as.data.frame(dam_vec)
dam_d  <- na.omit(dam_d)

write.table(dam_d, 'fg_Mur_damaging.txt',
            sep = '\t', col.names = FALSE, row.names = FALSE,
            quote = FALSE, na = 'NA')


GO <- read_delim("GO_Mur_damaging.txt", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)

colnames(GO) <- gsub(" ", "_", colnames(GO))

GO_f <- filter(GO, Fold_Enrichment >=2 & FDR< 0.05)

GO_f <- GO_f[order(GO_f$Fold_Enrichment, decreasing = TRUE),]  

write.table(GO_f, 'res_Mur_damaging.txt',
            sep = '\t', col.names = TRUE, row.names = FALSE,
            quote = FALSE, na = 'NA')