# EXTRACT ORTHOPAIRS FROM SYNGAP OUTPUT

###  C. orientalis & C. rubella

# Load libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# Load SynGAP output file
Co_Cr <- read_delim("syngap_res.txt", 
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE)

Co_Cr$Co <- substr( Co_Cr$X1, 1, 2)
Co_Cr$Cr <- substr( Co_Cr$X2, 1, 2)

# filter output files
data_filt <- filter(Co_Cr, Cr != "Cr" & Co != "Co" )

# select unique gene IDs
CO_unique <- as.data.frame(unique(data_filt$X1)) 
colnames(CO_unique) <- "X1"
CR_unique <-  as.data.frame(unique(data_filt$X2)) 
colnames(CR_unique) <- "X2"

# count occurrences of gene iDs
CO_count <- CO_unique %>%
  left_join(data_filt %>% group_by(X1) %>% summarise(count = n()), by = "X1") 

CR_count <- CR_unique %>%
  left_join(data_filt %>% group_by(X2) %>% summarise(count = n()), by = "X2") 

# extract orthopairs
CO_unique <- filter(CO_count, count==1) 
CR_unique <- filter(CR_count, count==1) 

CO_not_unique <- filter(CO_count, count!=1) 
CR__not_unique <- filter(CR_count, count!=1) 

CO_uni <- as.vector(CO_unique$X1)
CR_uni <- as.vector(CR_unique$X2)

orthopairs <- data_filt %>%
  filter(X1 %in% CO_uni & X2 %in% CR_uni)

orthopairs <- orthopairs[,c(1,2)]
length(unique(orthopairs$X1))
length(unique(orthopairs$X2))

#write.table(orthopairs, "Co_Cr_orthopairs.tsv", sep="\t", quote=FALSE, row.names = FALSE)


### Co, Cr with A. thaliana

Co <- read_delim("Co.Ath.final.genepair", 
                 delim = "\t", escape_double = FALSE, 
                 col_names = FALSE, trim_ws = TRUE)


Cr <- read_delim("Cr.Ath.final.genepair", 
                 delim = "\t", escape_double = FALSE, 
                 col_names = FALSE, trim_ws = TRUE)

orthopairs <- read_delim("../Co_cr_prot_syngap_chr_wise/orthopairs.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

orthopairs$Co <-  substr(orthopairs$Co, 1, nchar(orthopairs$Co) - 3)
orthopairs$Cr <-  substr(orthopairs$Cr, 1, nchar(orthopairs$Cr) - 3)

Co$Co <- substr( Co$X1, 1, 2)
Co$Ath <- substr( Co$X2, 1, 2)

Cr$Cr <- substr( Cr$X1, 1, 2)
Cr$Ath <- substr( Cr$X2, 1, 2)

Co_filt <- filter(Co, Co != "Co" & Ath != "At")
Cr_filt <- filter(Cr, Cr != "Cr" & Ath != "At")

# Co & Ath
CO_unique <- as.data.frame(unique(Co_filt$X1)) 
colnames(CO_unique) <- "X1"
CO_Ath_unique <-  as.data.frame(unique(Co_filt$X2))
colnames(CO_Ath_unique) <- "X2"

CO_count <- CO_unique %>%
  left_join(Co_filt %>% group_by(X1) %>% summarise(count = n()), by = "X1") 

CO_Ath_count <- CO_Ath_unique %>%
  left_join(Co_filt %>% group_by(X2) %>% summarise(count = n()), by = "X2") 

CO_unique <- filter(CO_count, count==1)
CO_Ath_unique <- filter(CO_Ath_count, count==1)

CO_uni <- as.vector(CO_unique$X1)
CO_Ath_uni <- as.vector(CO_Ath_unique$X2)

Co_orthopairs <- Co_filt %>%
  filter(X1 %in% CO_uni & X2 %in% CO_Ath_uni)

Co_orthopairs <- Co_orthopairs[,c(1,2)]
length(unique(Co_orthopairs$X1))
length(unique(Co_orthopairs$X2)) # 17,367

colnames(Co_orthopairs) <- c("Co", "Ath")

#write.table(Co_orthopairs, "Co_Ath_orthopairs.tsv", sep="\t", quote=FALSE, row.names = FALSE)


# Cr & Ath
CR_unique <- as.data.frame(unique(Cr_filt$X1)) 
colnames(CR_unique) <- "X1"
CR_Ath_unique <-  as.data.frame(unique(Cr_filt$X2))
colnames(CR_Ath_unique) <- "X2"

CR_count <- CR_unique %>%
  left_join(Cr_filt %>% group_by(X1) %>% summarise(count = n()), by = "X1") 

CR_Ath_count <- CR_Ath_unique %>%
  left_join(Cr_filt %>% group_by(X2) %>% summarise(count = n()), by = "X2") 

CR_unique <- filter(CR_count, count==1)
CR_Ath_unique <- filter(CR_Ath_count, count==1)

CR_uni <- as.vector(CR_unique$X1)
CR_Ath_uni <- as.vector(CR_Ath_unique$X2)

Cr_orthopairs <- Cr_filt %>%
  filter(X1 %in% CR_uni & X2 %in% CR_Ath_uni)

Cr_orthopairs <- Cr_orthopairs[,c(1,2)]
length(unique(Cr_orthopairs$X1))
length(unique(Cr_orthopairs$X2)) #17,435

colnames(Cr_orthopairs) <- c("Cr", "Ath") 

#write.table(Cr_orthopairs, "Cr_Ath_orthopairs.tsv", sep="\t", quote=FALSE, row.names = FALSE)