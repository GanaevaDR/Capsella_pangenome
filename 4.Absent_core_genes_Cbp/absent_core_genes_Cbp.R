# Analysis of PacBio data in Capsella bursa-pastoris

library(readr)
library(dplyr)

rm(list=ls())
folder_path <- "C:/Users/Dasha/Desktop/masters/Cbp/breadth"
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

merged_df <- data.frame()

for (file in file_list) {
  df <- read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  
  temp_df <- df %>%
    select(4, 14) %>%  
    rename(key = X4, value = X14)  
  if (nrow(merged_df) == 0) {
    merged_df <- temp_df
  } else {
    merged_df <- full_join(merged_df, temp_df, by = "key", suffix = c("", paste0("_", basename(file))))
  }
}

colnames(merged_df)[-1] <- gsub("txt", "", list.files(path = folder_path, pattern = "*.txt"))
names(merged_df)[names(merged_df) == 'key'] <- 'gene'

colnames(merged_df) <- substr(colnames(merged_df), 1, nchar(colnames(merged_df))-1) 

Cbp_breadth <- merged_df
Cbp_breadth <- as.data.frame(Cbp_breadth)
rownames(Cbp_breadth) <- Cbp_breadth$gen
Cbp_breadth$gene <- NULL

write.table(Cbp_breadth, 'Cbp_mean_breadth.txt',
            sep = '\t', col.names = TRUE, row.names = TRUE,
            quote = FALSE, na = 'NA')

b_vec <- unlist(Cbp_breadth[, ])

Cbp_breadth$mean <- rowMeans(Cbp_breadth, na.rm = TRUE)

length(which(Cbp_breadth$mean == 0)) #316
length(which(Cbp_breadth$mean ==1 )) #~50,000

#distribution plots
library(ggplot2)
library(gridExtra)

# subgenome-wise
Cbp_breadth$mean <- NULL

breadth_O <- Cbp_breadth[1:27095,]
breadth_R <- Cbp_breadth[27096:54186,]

o_vec <- unlist(breadth_O[,])
r_vec <- unlist(breadth_R[,])

r_vec_df <- as.data.frame(r_vec)
o_vec_df <- as.data.frame(o_vec)

r_all <- ggplot(r_vec_df, aes(x = r_vec)) +
  geom_histogram(bins = 25, 
                 fill = "blue", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "5 populations of C. bursa-pastoris (R)", 
       x = "Average gene breadth", 
       y = "Number of genes") +
  theme_minimal()


o_all <- ggplot(o_vec_df, aes(x = o_vec)) +
  geom_histogram(bins = 25, 
                 fill = "firebrick", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "5 populations of C. bursa-pastoris (O)", 
       x = "average gene breadth", 
       y = "Number of genes") +
  theme_minimal()

grid.arrange(r_all, o_all, nrow = 1, ncol = 2)


# work with processed files

library(readr)
library(dplyr)

P_bre <- read.table("Cbp_mean_breadth.txt", header = TRUE, row.names = 1,
                    sep="\t", na.strings = "NA", dec = ".",
                    strip.white = TRUE)

conserved <- read_delim("conserved_genes_filt.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

conserved <- conserved[,1:2]
cons_vec <- c(conserved$Co, conserved$Cr) 

v1 <- rownames(P_bre[P_bre$res_Eng==0 ,]) 
v2 <- rownames(P_bre[P_bre$res_KBG==0 ,]) 
v3 <- rownames(P_bre[P_bre$res_le3==0 ,]) 
v4 <- rownames(P_bre[P_bre$res_MSK==0 ,])
v5 <- rownames(P_bre[P_bre$res_Mur==0 ,])

a <- intersect(v1, v2)
a1 <- intersect(a,v3)
a2 <- intersect(a1, v4)
a3 <- intersect(a2, v5) #316

abs_Cbp <- intersect(a3, cons_vec)
abs_df <- as.data.frame(abs_Cbp)

#write.table(abs_df, 'Cbp_PacBio_bre_1_conserved.txt',
#                   sep = '\t', col.names = FALSE, row.names = FALSE,
#            quote = FALSE, na = 'NA')

library(VennDiagram)
library(grid)

# Combine the sets into a list
venn_list <- list(Set1 = v1, Set2 = v2, Set3 = v3, Set4 = v4, Set5 = v5)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Eng", "KBG", "le3", "MSK", "Mur"),
  filename = NULL, 
  output = TRUE,
  lwd = 2,
  lty = 'solid',
  fill = c("red", "blue", "green", "yellow", "purple"),
  alpha = 0.4,
  label.col = "black",
  cex = 1,
  fontfamily = "serif",
  cat.cex = 1,
  cat.fontfamily = "serif",
  
  margin = 0.1, 
  height = 10,
  width = 10
)

# Plot the Venn diagram
grid.draw(venn.plot)
grid.text("Genes with mean breadth 0", x = 0.5, y = 0.95, gp = gpar(fontsize = 10, fontface = "bold"))

length(intersect(a3, cons_vec)) #14 from 316 are conserved in parent species

# Analyze conserved genes with breadth 1

v1 <- rownames(P_bre[P_bre$res_Eng==1 ,]) 
v2 <- rownames(P_bre[P_bre$res_KBG==1 ,]) 
v3 <- rownames(P_bre[P_bre$res_le3==1 ,]) 
v4 <- rownames(P_bre[P_bre$res_MSK==1 ,]) 
v5 <- rownames(P_bre[P_bre$res_Mur==1 ,])

a <- intersect(v1, v2)
a1 <- intersect(a,v3)
a2 <- intersect(a1, v4)
a3 <- intersect(a2, v5) #49k 

length(intersect(a3, cons_vec)) #39k belong to conserved parent genes

cons_Cbp <- intersect(a3, cons_vec)
cons_df <- as.data.frame(cons_Cbp)

#write.table(cons_df, 'Cbp_PacBio_bre_1_conserved.txt',
 #                   sep = '\t', col.names = FALSE, row.names = FALSE,
  #            quote = FALSE, na = 'NA')

library(VennDiagram)
library(grid)

venn_list <- list(Set1 = v1, Set2 = v2, Set3 = v3, Set4 = v4, Set5 = v5)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Eng", "KBG", "le3", "MSK", "Mur"),
  filename = NULL, 
  output = TRUE,
  lwd = 2,
  lty = 'solid',
  fill = c("red", "blue", "orange", "green", "purple"),
  alpha = 0.4,
  label.col = "black",
  cex = 1,
  fontfamily = "serif",
  cat.cex = 1,
  cat.fontfamily = "serif",
  
  margin = 0.1,
  height = 10, 
  width = 10
)

# Plot the Venn diagram
grid.draw(venn.plot)
grid.text("Genes with average gene breadth 1", x = 0.5, y = 0.95, gp = gpar(fontsize = 10, fontface = "bold"))


#=================================================================================
#=================================================================================

# # #              1/0 fraction in orthopairs
P_bre <- read.table("Cbp_mean_breadth.txt", header = TRUE, row.names = 1,
                    sep="\t", na.strings = "NA", dec = ".",
                    strip.white = TRUE)

v1 <- rownames(P_bre[P_bre$res_Eng==0 ,]) #1401
v2 <- rownames(P_bre[P_bre$res_KBG==0 ,]) #1514
v3 <- rownames(P_bre[P_bre$res_le3==0 ,]) #1336
v4 <- rownames(P_bre[P_bre$res_MSK==0 ,]) #1497
v5 <- rownames(P_bre[P_bre$res_Mur==0 ,]) #1376

v1_df <- as.data.frame(v1)
v2_df <- as.data.frame(v2)
v3_df <- as.data.frame(v3)
v4_df <- as.data.frame(v4)
v5_df <- as.data.frame(v5)

a <- intersect(v1, v2)
a1 <- intersect(a,v3)
a2 <- intersect(a1, v4)
a3 <- intersect(a2, v5) 

abs_Cbp <- a3 
abs_Cbp <- intersect(a3, cons_vec)

# line 1
v1_O <- as.data.frame(v1_df[1:562,])
colnames(v1_O) <- "Co"
v1_O$v1_O_abs <- 1

v1_R <- as.data.frame(v1_df[563:1401,])
colnames(v1_R) <- "Cr"
v1_R$v1_R_abs <- 1

# line 2
v2_O <- as.data.frame(v2_df[1:262,])
colnames(v2_O) <- "Co"
v2_O$v2_O_abs <- 1

v2_R <- as.data.frame(v2_df[263:1514,])
colnames(v2_R) <- "Cr"
v2_R$v2_R_abs <- 1

# line 3
v3_O <- as.data.frame(v3_df[1:410,])
colnames(v3_O) <- "Co"
v3_O$v3_O_abs <- 1

v3_R <- as.data.frame(v3_df[411:1336,])
colnames(v3_R) <- "Cr"
v3_R$v3_R_abs <- 1

# line 4
v4_O <- as.data.frame(v4_df[1:562,])
colnames(v4_O) <- "Co"
v4_O$v4_O_abs <- 1

v4_R <- as.data.frame(v4_df[563:1497,])
colnames(v4_R) <- "Cr"
v4_R$v4_R_abs <- 1

# line 5
v5_O <- as.data.frame(v5_df[1:465,])
colnames(v5_O) <- "Co"
v5_O$v5_O_abs <- 1

v5_R <- as.data.frame(v5_df[466:1376,])
colnames(v5_R) <- "Cr"
v5_R$v5_R_abs <- 1

v1_abs_O <- merge(conserved, v1_O, by="Co", all.x = TRUE)
v1_abs_OR <- merge(v1_abs_O, v1_R, by="Cr", all.x = TRUE)
v1_abs_OR[is.na(v1_abs_OR)] <- 0

v2_abs_O <- merge(conserved, v2_O, by="Co", all.x = TRUE)
v2_abs_OR <- merge(v2_abs_O, v2_R, by="Cr", all.x = TRUE)
v2_abs_OR[is.na(v2_abs_OR)] <- 0
v2_abs_OR$Cr <- NULL

v3_abs_O <- merge(conserved, v3_O, by="Co", all.x = TRUE)
v3_abs_OR <- merge(v3_abs_O, v3_R, by="Cr", all.x = TRUE)
v3_abs_OR[is.na(v3_abs_OR)] <- 0
v3_abs_OR$Cr <- NULL

v4_abs_O <- merge(conserved, v4_O, by="Co", all.x = TRUE)
v4_abs_OR <- merge(v4_abs_O, v4_R, by="Cr", all.x = TRUE)
v4_abs_OR[is.na(v4_abs_OR)] <- 0
v4_abs_OR$Cr <- NULL

v5_abs_O <- merge(conserved, v5_O, by="Co", all.x = TRUE)
v5_abs_OR <- merge(v5_abs_O, v5_R, by="Cr", all.x = TRUE)
v5_abs_OR[is.na(v5_abs_OR)] <- 0
v5_abs_OR$Cr <- NULL

m1 <- merge(v1_abs_OR, v2_abs_OR, by="Co")
m2 <- merge(m1,v3_abs_OR, by="Co" )
m3 <- merge(m2,v4_abs_OR, by="Co" )
m4 <- merge(m3,v5_abs_OR, by="Co" )

# filtration
f1 <- filter(m4, ((v1_O_abs == 1 & v1_R_abs != 1) | (v1_O_abs != 1 & v1_R_abs == 1)) | ((v2_O_abs == 1 & v2_R_abs != 1) | (v2_O_abs != 1 & v2_R_abs == 1)) |((v3_O_abs == 1 & v3_R_abs != 1) | (v3_O_abs != 1 & v3_R_abs == 1)) |((v4_O_abs == 1 & v4_R_abs != 1) | (v4_O_abs != 1 & v4_R_abs == 1)) |((v5_O_abs == 1 & v5_R_abs != 1) | (v5_O_abs != 1 & v5_R_abs == 1)))

fin_abs <- f1

fin_abs$sum_O <- rowSums(fin_abs[,c(3,5,7,9,11)])
fin_abs$sum_R <- rowSums(fin_abs[,c(4,6,8,10,12)])
fin_abs$sum_all <- rowSums(fin_abs[,c(3:12)])


#write.table(fin_abs, 'C:/Users/Dasha/Desktop/masters/IITP/pop_data_qc/Cbp/PacBio/pb_MQ_60/01_fraction_genes_Cbp.txt',
 #          sep = '\t', col.names = TRUE, row.names = FALSE,
  #        quote = FALSE, na = 'NA')

#=================================================================================
#=================================================================================

# # #                               DRAW PLOTS

rm(list=ls())
library(readr)
library(dplyr)

P_bre <- read.table("Cbp_mean_breadth.txt", header = TRUE, row.names = 1,
                    sep="\t", na.strings = "NA", dec = ".",
                    strip.white = TRUE)
P_cov$mean <- NULL
P_vec <- unlist(P_cov[,])
P_vec_df <- as.data.frame(P_vec)

library(ggplot2)
library(gridExtra)

ggplot(P_vec_df, aes(x = P_vec)) +
  geom_histogram(bins = 100, 
                 fill = "blue", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "", 
       x = "Average read breadth", 
       y = "Number of genes") +
  theme_minimal()

coverage_O <- P_cov[1:27095,]
coverage_R <- P_cov[27096:54186,]

o_vec <- unlist(coverage_O[,])
r_vec <- unlist(coverage_R[,])

r_vec_df <- as.data.frame(r_vec)
o_vec_df <- as.data.frame(o_vec)

r_all <- ggplot(r_vec_df, aes(x = r_vec)) +
  geom_histogram(bins = 50, 
                 fill = "blue", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "5 samples of Capsella bursa-pastoris (R)", 
       x = "Average read breadth", 
       y = "Number of genes") +
  scale_x_continuous(limits = c(-10,100))+
  theme_minimal()


o_all <- ggplot(o_vec_df, aes(x = o_vec)) +
  geom_histogram(bins = 50, 
                 fill = "firebrick", 
                 color = "black", 
                 alpha = 0.7) +  
  labs(title = "5 samples of Capsella bursa-pastoris (O)", 
       x = "Average read breadth", 
       y = "Number of genes") +
  scale_x_continuous(limits = c(-10,100))+
  theme_minimal()

grid.arrange(r_all, o_all, nrow = 1, ncol = 2)