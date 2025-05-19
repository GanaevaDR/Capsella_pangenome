# Load libraries
library(readr)
library(ggplot2)
library(ggdendro)
library(DESeq2)
library(PCAtools)
library(matrixStats)
library(ggplot2)
library(reshape2) 

RNA_seq_all_readcounts <- read_delim("RNA_seq_all_readcounts.csv", 
                                     delim = ";", escape_double = FALSE, trim_ws = TRUE)

RNA_seq_all_readcounts <- as.data.frame(RNA_seq_all_readcounts)
rownames(RNA_seq_all_readcounts) <- RNA_seq_all_readcounts$gene_id
RNA_seq_all_readcounts$gene_id <- NULL

normal_expr <- RNA_seq_all_readcounts[,1:63]

correlation <- cor(normal_expr)
correlation

distance <- as.dist(1 - correlation)
tree <- hclust(distance)

#draw the tree of hierarchical clustering
ggdendrogram(tree, rotate = TRUE, theme_dendro = FALSE) +
  labs(x = '', y = "1 - corr") + 
  theme_classic() + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

expr.matrix <- as.matrix(normal_expr)
expr.design <- data.frame(row.names =
                            colnames(expr.matrix),
                          condition = rep("a", ncol(normal_expr)))

# normalization using DESeq2
dds <- DESeqDataSetFromMatrix(countData = expr.matrix,
                              colData = expr.design,
                              design = ~ 1)

dds <- DESeq(dds)
res <- results(dds)
dds <- estimateSizeFactors(dds)
norm.expr <- counts(dds, normalized = TRUE)

#write.table(norm.expr, "Norm_expression_controls.txt",
 #           sep = "\t", col.names = TRUE, row.names = TRUE,
  #          quote = FALSE, na = "NA")

# log-scaling
norm.expr <- norm.expr+0.01
log.norm.expr <- log10(norm.expr)

#write.table(log.norm.expr, "Log_norm_expression_controls.txt",
 #          sep = "\t", col.names = TRUE, row.names = TRUE,
  #        quote = FALSE, na = "NA")

log_expr <- as.data.frame(log.norm.expr)
log_expr$gene_id <- rownames(log.norm.expr)
rownames(log_expr) <- NULL

# Group by gene fractions
# conserved genes
conserved <- read_csv("Cbp_covered_among_conserved_parent.txt", 
                                               col_names = FALSE)
colnames(conserved) <- "gene_id"

# absent genes
absent_all <- read_csv("Cbp_absent_all.txt", 
                       col_names = FALSE)
colnames(absent_all) <- "gene_id"

#1/0 fraction
X1_0_fraction <- read.table("10_fraction_among_conserved_parent.txt", 
                             header = FALSE, sep = "\t")
X1_0_fraction <- X1_0_fraction[,1:2]
X1_0_frac <- as.data.frame(c(X1_0_fraction$V1,X1_0_fraction$V2))
colnames(X1_0_frac) <- "gene_id"

# mutations (for each population)
mutations_Eng <- read_csv("Eng_damaged_exons.txt", 
                          col_names = FALSE)
mutations_Eng <- as.data.frame(unique(mutations_Eng$X1)) # 4068 ген
colnames(mutations_Eng) <- "gene_id"

mutations_KBG <- read_csv("KBG_damaged_exons.txt", 
                          col_names = FALSE)
mutations_KBG <- as.data.frame(unique(mutations_KBG$X1)) # 4425
colnames(mutations_KBG) <- "gene_id"

mutations_MSK <- read_csv("MSK_damaged_exons.txt", 
                          col_names = FALSE)
mutations_MSK <- as.data.frame(unique(mutations_MSK$X1)) # 3902
colnames(mutations_MSK) <- "gene_id"

# Group by populations
Eng <- log_expr[,c(64, 1:18)]
KBG <- log_expr[,c(64,19:35)]
MSK <- log_expr[,c(64,36:63)]

# merge with expression data population-wise
Eng_conserved <- merge(conserved, Eng, by="gene_id")
Eng_absent <- merge(absent, Eng, by="gene_id")
Eng_absent_all <- merge(absent_all, Eng, by="gene_id")
Eng_1_0 <- merge(X1_0_frac, Eng, by="gene_id")
Eng_mut <- merge(mutations_Eng, Eng, by="gene_id")

KBG_conserved <- merge(conserved, KBG, by="gene_id")
KBG_absent <- merge(absent, KBG, by="gene_id")
KBG_absent_all <- merge(absent_all, KBG, by="gene_id")
KBG_1_0 <- merge(X1_0_frac, KBG, by="gene_id")
KBG_mut <- merge(mutations_KBG, KBG, by="gene_id")

MSK_conserved <- merge(conserved, MSK, by="gene_id")
MSK_absent <- merge(absent, MSK, by="gene_id")
MSK_absent_all <- merge(absent_all, MSK, by="gene_id")
MSK_1_0 <- merge(X1_0_frac, MSK, by="gene_id")
MSK_mut <- merge(mutations_MSK, MSK, by="gene_id")

# calculate mean, median, max, sum for all groups of genes and populations
#Eng
Eng_conserved$mean <- rowMeans(Eng_conserved[2:19],  na.rm = TRUE)
Eng_conserved$median <- apply(Eng_conserved[2:19], 1, median)
Eng_conserved$max <- apply(Eng_conserved[2:19], 1, max)
Eng_conserved$sum <- apply(Eng_conserved[2:19], 1, sum)

Eng_absent_all$mean <- rowMeans(Eng_absent_all[2:19],  na.rm = TRUE)
Eng_absent_all$median <- apply(Eng_absent_all[2:19], 1, median)
Eng_absent_all$max <- apply(Eng_absent_all[2:19], 1, max)
Eng_absent_all$sum <- apply(Eng_absent_all[2:19], 1, sum)

Eng_1_0$mean <- rowMeans(Eng_1_0[2:19],  na.rm = TRUE)
Eng_1_0$median <- apply(Eng_1_0[2:19], 1, median)
Eng_1_0$max <- apply(Eng_1_0[2:19], 1, max)
Eng_1_0$sum <- apply(Eng_1_0[2:19], 1, sum)

Eng_mut$mean <- rowMeans(Eng_mut[2:19],  na.rm = TRUE)
Eng_mut$median <- apply(Eng_mut[2:19], 1, median)
Eng_mut$max <- apply(Eng_mut[2:19], 1, max)
Eng_mut$sum <- apply(Eng_mut[2:19], 1, sum)

# KBG
KBG_conserved$mean <- rowMeans(KBG_conserved[2:18],  na.rm = TRUE)
KBG_conserved$median <- apply(KBG_conserved[2:18], 1, median)
KBG_conserved$max <- apply(KBG_conserved[2:18], 1, max)
KBG_conserved$sum <- apply(KBG_conserved[2:18], 1, sum)

KBG_absent_all$mean <- rowMeans(KBG_absent_all[2:18],  na.rm = TRUE)
KBG_absent_all$median <- apply(KBG_absent_all[2:18], 1, median)
KBG_absent_all$max <- apply(KBG_absent_all[2:18], 1, max)
KBG_absent_all$sum <- apply(KBG_absent_all[2:18], 1, sum)

KBG_1_0$mean <- rowMeans(KBG_1_0[2:18],  na.rm = TRUE)
KBG_1_0$median <- apply(KBG_1_0[2:18], 1, median)
KBG_1_0$max <- apply(KBG_1_0[2:18], 1, max)
KBG_1_0$sum <- apply(KBG_1_0[2:18], 1, sum)

KBG_mut$mean <- rowMeans(KBG_mut[2:18],  na.rm = TRUE)
KBG_mut$median <- apply(KBG_mut[2:18], 1, median)
KBG_mut$max <- apply(KBG_mut[2:18], 1, max)
KBG_mut$sum <- apply(KBG_mut[2:18], 1, sum)

# MSK
MSK_conserved$mean <- rowMeans(MSK_conserved[2:29],  na.rm = TRUE)
MSK_conserved$median <- apply(MSK_conserved[2:29], 1, median)
MSK_conserved$max <- apply(MSK_conserved[2:29], 1, max)
MSK_conserved$sum <- apply(MSK_conserved[2:29], 1, sum)

MSK_absent_all$mean <- rowMeans(MSK_absent_all[2:29],  na.rm = TRUE)
MSK_absent_all$median <- apply(MSK_absent_all[2:29], 1, median)
MSK_absent_all$max <- apply(MSK_absent_all[2:29], 1, max)
MSK_absent_all$sum <- apply(MSK_absent_all[2:29], 1, sum)

MSK_1_0$mean <- rowMeans(MSK_1_0[2:29],  na.rm = TRUE)
MSK_1_0$median <- apply(MSK_1_0[2:29], 1, median)
MSK_1_0$max <- apply(MSK_1_0[2:29], 1, max)
MSK_1_0$sum <- apply(MSK_1_0[2:29], 1, sum)

MSK_mut$mean <- rowMeans(MSK_mut[2:29],  na.rm = TRUE)
MSK_mut$median <- apply(MSK_mut[2:29], 1, median)
MSK_mut$max <- apply(MSK_mut[2:29], 1, max)
MSK_mut$sum <- apply(MSK_mut[2:29], 1, sum)

#Eng
# Select columns for visualization
df <- Eng_conserved[,20:23]
df <- Eng_absent_all[,20:23]
df <- Eng_1_0[,20:23]
df <- Eng_mut[,20:23]

length(which(Eng_conserved$mean == -2)) #458
length(which(Eng_absent_all$mean == -2)) #183
length(which(Eng_1_0$mean == -2)) #85
length(which(Eng_mut$mean == -2)) #200

# Reshape the dataframe
df_long <- melt(df)

colors <- c("mean" = "#A8DADC", "median" = "#457B9D", 
            "max" = "#F8BBD0", "sum" = "#B2F2BB")

# Create distribution plots 
ggplot(df_long, aes(x = value, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.7, position = 'identity', 
                 color = "black", 
                 linewidth = 0.5) + 
  facet_wrap(~ variable, ncol = 4, scales = 'free') +  
  scale_fill_manual(values = colors) + 
  labs(title = "Distribution of expression level of C. bursa-pastoris genes (ENG)",
       x = "Normalized expression level, log10",
       y = "Number of genes") +
  theme_minimal()

#===============================================================================

df <- KBG_conserved[,19:22]
df <- KBG_absent_all[,19:22]
df <- KBG_1_0[,19:22]
df <- KBG_mut[,19:22]

length(which(KBG_conserved$mean == -2)) #625
length(which(KBG_absent_all$mean == -2)) #223
length(which(KBG_1_0$mean == -2)) #127
length(which(KBG_mut$mean == -2)) #270

df_long <- melt(df)

colors <- c("mean" = "#A8DADC", "median" = "#457B9D", 
            "max" = "#F8BBD0", "sum" = "#B2F2BB")

ggplot(df_long, aes(x = value, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.7, position = 'identity', 
                 color = "black", 
                 linewidth = 0.5) + 
  facet_wrap(~ variable, ncol = 4, scales = 'free') +  
  scale_fill_manual(values = colors) +  
  labs(title = "Distribution of expression level of C. bursa-pastoris genes (KBG)",
       x = "Normalized expression level, log10",
       y = "Number of genes") +
  theme_minimal()

#===============================================================================

df <- MSK_conserved[,30:33]
df <- MSK_absent_all[,30:33]
df <- MSK_1_0[,30:33]
df <- MSK_mut[,30:33]

length(which(MSK_conserved$mean == -2)) #516
length(which(MSK_absent_all$mean == -2)) #207
length(which(MSK_1_0$mean == -2)) #106
length(which(MSK_mut$mean == -2)) #223

df_long <- melt(df)

colors <- c("mean" = "#A8DADC", "median" = "#457B9D", 
            "max" = "#F8BBD0", "sum" = "#B2F2BB")

ggplot(df_long, aes(x = value, fill = variable)) +
  geom_histogram(bins = 30, alpha = 0.7, position = 'identity', 
                 color = "black", 
                 linewidth = 0.5) + 
  facet_wrap(~ variable, ncol = 4, scales = 'free') +  
  scale_fill_manual(values = colors) +  
  labs(title = "Distribution of expression level of C. bursa-pastoris genes (MSK)",
       x = "Normalized expression level, log10",
       y = "Number of genes") +
  theme_minimal()