library("ampir")
library(readr)
library(dplyr)

Co_cds <- read_faa("Co_cds.faa")
Cr_cds <- read_faa("Cr_cds.faa")

orthopairs <- read_delim("orthopairs.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

orthopairs$Co <-  substr(orthopairs$Co, 1, nchar(orthopairs$Co) - 3)
orthopairs$Cr <-  substr(orthopairs$Cr, 1, nchar(orthopairs$Cr) - 3)

colnames(orthopairs)=c("Co_gene_name","Cr_gene_name")

colnames(Co_cds)=c("Co_gene","Co_seq")
colnames(Cr_cds)=c("Cr_gene","Cr_seq")

Co_cds$Co_len <- nchar(Co_cds$Co_seq)
Cr_cds$Cr_len <- nchar(Cr_cds$Cr_seq)

Co_cds$name <- gsub( " .*$", "", Co_cds$Co_gene)
Co_cds$Co_gene_name <- sub("\\.[^.]*$", "", Co_cds$name)
Co_cds$name = NULL

Cr_cds$name <- gsub( " .*$", "", Cr_cds$Cr_gene)
Cr_cds$Cr_gene_name <- sub("\\.[^.]*$", "", Cr_cds$name)
Cr_cds$name = NULL

Co_Cr <- merge(orthopairs, Co_cds, by="Co_gene_name")
Co_Cr <- merge(Co_Cr, Cr_cds, by="Cr_gene_name")

sum(is.na(Co_Cr$Co_gene_name))

Co_Cr$diff <- Co_Cr$Co_len - Co_Cr$Cr_len

Co_Cr_diff <- filter(Co_Cr, diff !=0) 

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"Cr_gene_name"], "_", data[rowNum,"abs_diff"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Cr_seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

Cr_sel <- filter(Co_Cr, diff < 0)
Cr_sel <- filter(Co_Cr, Cr_len > Co_len)
Cr_sel$abs_diff <- abs(Cr_sel$diff)

writeFasta(Cr_sel, "C:/Users/Dasha/Desktop/masters/IITP/minimap/new_anno/Cr_larger.fasta")


writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"Co_gene_name"], "_", data[rowNum,"diff"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Co_seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

Co_sel <- filter(Co_Cr, diff > 0)
Co_sel <- filter(Co_Cr, Co_len > Cr_len)
writeFasta(Co_sel, "C:/Users/Dasha/Desktop/masters/IITP/minimap/new_anno/Co_larger.fasta")

Co_Cr_diff <- Co_Cr[, c(1,8,2,5)]
Co_Cr_diff <- filter(Co_Cr_diff, Co_len != Cr_len)
Co_Cr_diff$Cr_sh <- 0
Co_Cr_diff$Co_sh <- 0

Co_Cr_diff$Cr_sh <- ifelse(Co_Cr_diff[, 2] < Co_Cr_diff[, 4], 1, 0)
Co_Cr_diff$Co_sh <- ifelse(Co_Cr_diff[, 2] > Co_Cr_diff[, 4], 1, 0)

colnames(Co_Cr_diff) = c("Cr_gene",	"Cr_gene_length",	"Co_gene",	"Co_gene_length",	"Cr_sh"	,"Co_sh")
write.table(Co_Cr_diff, "orthopairs_diff.txt", sep = "\t", quote = FALSE,  row.names = FALSE)