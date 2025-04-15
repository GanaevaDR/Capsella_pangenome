# draw plots for CDS length difference distribution

library("ampir")
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

Co_cds <- read_faa("Co_cds.faa")
Cr_cds <- read_faa("Cr_cds.faa")

orthopairs <- read_delim("orthopairs.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

orthopairs$Co <-  substr(orthopairs$Co, 1, nchar(orthopairs$Co) - 3)
orthopairs$Cr <-  substr(orthopairs$Cr, 1, nchar(orthopairs$Cr) - 3)
colnames(orthopairs) <- c("Co_gene_name","Cr_gene_name")

colnames(Co_cds)=c("Co_gene","Co_seq")
colnames(Cr_cds)=c("Cr_gene","Cr_seq")

Co_cds$Co_len <- nchar(Co_cds$Co_seq)
Cr_cds$Cr_len <- nchar(Cr_cds$Cr_seq)

Co_merged <- merge(orthopairs, Co_cds, by="Co_gene_name")

ortho_merged <- merge(Co_merged, Cr_cds, by="Cr_gene_name")
ortho_merged$d_val <- ortho_merged$Co_len - ortho_merged$Cr_len

Co_larger <- filter(ortho_merged, Co_len >= Cr_len)
Cr_larger <- filter(ortho_merged, Cr_len >= Co_len)

Cr_larger$d_val <- Cr_larger$Cr_len - Cr_larger$Co_len

Co_larger_mod <- Co_larger
Cr_larger_mod <- Cr_larger

Co_larger_mod[Co_larger_mod$d_val > 50,]$d_val <- 50
Cr_larger_mod[Cr_larger_mod$d_val > 50,]$d_val <- 50


p1 <- ggplot(Co_larger_mod, aes(x = d_val)) +
  geom_histogram(bins = 50, 
                 fill = "blue", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "", 
       x = "Difference of CDS length (longer C. orientalis CDSs)", 
       y = "Number of CDSs") +
  theme_minimal()

p2 <- ggplot(Cr_larger_mod, aes(x = d_val)) +
  geom_histogram(bins = 50, 
                 fill = "lightblue", 
                 color = "black", 
                 alpha = 0.7) + 
  labs(title = "", 
       x = "Difference of CDS length (longer C. rubella CDSs)", 
       y = "Number of CDSs") +
  theme_minimal()

grid.arrange(p1, p2, nrow = 1, ncol = 2)