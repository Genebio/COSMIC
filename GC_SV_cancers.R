setwd("C:/Users/user/Downloads/COSMIC")
library(dplyr)
library(openxlsx)

#SV updated different GC types: diffuse, intestinal, mixed
clinical_data <- read.xlsx("Updated_GC_SV_wang2014.xlsx", sheet=1)
str(clinical_data)
SV_data <- read.xlsx("Updated_GC_SV_wang2014.xlsx", sheet=2)
str(SV_data)
clinical_data_matched <- data.frame(clinical_data[match(SV_data$Sample, clinical_data$Sample),2:3], row.names=NULL)
chrom <- paste0("chr", SV_data$Chr)
GC_SV_merged_df <- data.frame(chrom, Start=SV_data$Start, End=SV_data$End, Group=clinical_data_matched$GC_type, SV_type=SV_data$SV, Strand=SV_data$Strand)
write.table(GC_SV_merged_df, file="GC_SV_merged_df.bed", sep="\t", quote = F, row.names = F, col.names = F)

GC_SV_intestinal <- GC_SV_merged_df %>% filter(Group=="intestinal")
write.table(GC_SV_intestinal, file="GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
GC_SV_diffuse <- GC_SV_merged_df %>% filter(Group=="diffuse")
write.table(GC_SV_diffuse, file="GC_SV_diffuse.bed", sep="\t", quote = F, row.names = F, col.names = F)
GC_SV_mixed <- GC_SV_merged_df %>% filter(Group=="mixed")
write.table(GC_SV_mixed, file="GC_SV_mixed.bed", sep="\t", quote = F, row.names = F, col.names = F)

GC_merged_types_gene_sizes_intersect <- readr::read_tsv("GC_merged_types_gene_sizes_intersect.bed", col_names = F) %>% distinct() %>% filter(X8!=0)
#table(GC_merged_types_gene_sizes_intersect$X8)
nrow(GC_merged_types_gene_sizes_intersect)
gene_length=GC_merged_types_gene_sizes_intersect$X3-GC_merged_types_gene_sizes_intersect$X2
gene <- GC_merged_types_gene_sizes_intersect$X4
GC_group <- GC_merged_types_gene_sizes_intersect$X7
GC_SV_sizes_df <- data.frame(gene, gene_length, GC_group)
View(GC_SV_sizes_df)
GC_SV_sizes_df <- GC_SV_sizes_df %>% group_by(gene, GC_group) %>% summarise(gene_length=sum(gene_length))
GC_SV_sizes_df$gene_length <- log(GC_SV_sizes_df$gene_length)
library("RColorBrewer")
pdf("GC_SV_sizes_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~GC_group, data = GC_SV_sizes_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
summary(aov(gene_length~GC_group, data = GC_SV_sizes_df))
table(GC_SV_sizes_df$GC_group)

table(GC_SV_merged_df$SV_type)
#CTX
CTX_GC_SV_merged_df <- GC_SV_merged_df %>% filter(SV_type=="CTX")
CTX_GC_SV_intestinal <- CTX_GC_SV_merged_df %>% filter(Group=="intestinal")
write.table(CTX_GC_SV_intestinal, file="CTX_GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
CTX_GC_SV_diffuse <- CTX_GC_SV_merged_df %>% filter(Group=="diffuse")
write.table(CTX_GC_SV_diffuse, file="CTX_GC_SV_diffuse.bed", sep="\t", quote = F, row.names = F, col.names = F)
CTX_GC_SV_mixed <- CTX_GC_SV_merged_df %>% filter(Group=="mixed")
write.table(CTX_GC_SV_mixed, file="CTX_GC_SV_mixed.bed", sep="\t", quote = F, row.names = F, col.names = F)

CTX_GC_merged_types_gene_sizes_intersect <- readr::read_tsv("CTX_GC_merged_types_gene_sizes_intersect.bed", col_names = F) %>% distinct() %>% filter(X8!=0)
#table(GC_merged_types_gene_sizes_intersect$X8)
nrow(CTX_GC_merged_types_gene_sizes_intersect)
gene_length=CTX_GC_merged_types_gene_sizes_intersect$X3-CTX_GC_merged_types_gene_sizes_intersect$X2
gene <- CTX_GC_merged_types_gene_sizes_intersect$X4
GC_group <- CTX_GC_merged_types_gene_sizes_intersect$X7
GC_SV_sizes_df <- data.frame(gene, gene_length, GC_group)
#View(GC_SV_sizes_df)
GC_SV_sizes_df <- GC_SV_sizes_df %>% group_by(gene, GC_group) %>% summarise(gene_length=sum(gene_length))
GC_SV_sizes_df$gene_length <- log(GC_SV_sizes_df$gene_length)
library("RColorBrewer")
pdf("CTX_GC_SV_sizes_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~GC_group, data = GC_SV_sizes_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
summary(aov(gene_length~GC_group, data = GC_SV_sizes_df))
CTX_GC_stat <- data.frame(GC_SV_sizes_df) %>% wilcox_test(gene_length~GC_group)
write.xlsx(CTX_GC_stat, file="CTX_GC_stat.xlsx")

#DEL
DEL_GC_SV_merged_df <- GC_SV_merged_df %>% filter(SV_type=="DEL")
DEL_GC_SV_intestinal <- DEL_GC_SV_merged_df %>% filter(Group=="intestinal")
write.table(DEL_GC_SV_intestinal, file="DEL_GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
DEL_GC_SV_diffuse <- DEL_GC_SV_merged_df %>% filter(Group=="diffuse")
write.table(DEL_GC_SV_diffuse, file="DEL_GC_SV_diffuse.bed", sep="\t", quote = F, row.names = F, col.names = F)
DEL_GC_SV_mixed <- DEL_GC_SV_merged_df %>% filter(Group=="mixed")
write.table(DEL_GC_SV_mixed, file="DEL_GC_SV_mixed.bed", sep="\t", quote = F, row.names = F, col.names = F)

DEL_GC_merged_types_gene_sizes_intersect <- readr::read_tsv("DEL_GC_merged_types_gene_sizes_intersect.bed", col_names = F) %>% distinct() %>% filter(X8!=0)
#table(GC_merged_types_gene_sizes_intersect$X8)
nrow(DEL_GC_merged_types_gene_sizes_intersect)
gene_length=DEL_GC_merged_types_gene_sizes_intersect$X3-DEL_GC_merged_types_gene_sizes_intersect$X2
gene <- DEL_GC_merged_types_gene_sizes_intersect$X4
GC_group <- DEL_GC_merged_types_gene_sizes_intersect$X7
GC_SV_sizes_df <- data.frame(gene, gene_length, GC_group)
#View(GC_SV_sizes_df)
GC_SV_sizes_df <- GC_SV_sizes_df %>% group_by(gene, GC_group) %>% summarise(gene_length=sum(gene_length))
GC_SV_sizes_df$gene_length <- log(GC_SV_sizes_df$gene_length)
library("RColorBrewer")
pdf("DEL_GC_SV_sizes_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~GC_group, data = GC_SV_sizes_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
summary(aov(gene_length~GC_group, data = GC_SV_sizes_df))
table(GC_SV_sizes_df$GC_group)
DEL_GC_stat <- data.frame(GC_SV_sizes_df) %>% wilcox_test(gene_length~GC_group)
write.xlsx(DEL_GC_stat, file="DEL_GC_stat.xlsx")

#ITX
ITX_GC_SV_merged_df <- GC_SV_merged_df %>% filter(SV_type=="ITX")
ITX_GC_SV_intestinal <- ITX_GC_SV_merged_df %>% filter(Group=="intestinal")
write.table(ITX_GC_SV_intestinal, file="ITX_GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
ITX_GC_SV_diffuse <- ITX_GC_SV_merged_df %>% filter(Group=="diffuse")
write.table(ITX_GC_SV_diffuse, file="ITX_GC_SV_diffuse.bed", sep="\t", quote = F, row.names = F, col.names = F)
ITX_GC_SV_mixed <- ITX_GC_SV_merged_df %>% filter(Group=="mixed")
write.table(ITX_GC_SV_mixed, file="ITX_GC_SV_mixed.bed", sep="\t", quote = F, row.names = F, col.names = F)

ITX_GC_merged_types_gene_sizes_intersect <- readr::read_tsv("ITX_GC_merged_types_gene_sizes_intersect.bed", col_names = F) %>% distinct() %>% filter(X8!=0)
#table(GC_merged_types_gene_sizes_intersect$X8)
nrow(ITX_GC_merged_types_gene_sizes_intersect)
gene_length=ITX_GC_merged_types_gene_sizes_intersect$X3-ITX_GC_merged_types_gene_sizes_intersect$X2
gene <- ITX_GC_merged_types_gene_sizes_intersect$X4
GC_group <- ITX_GC_merged_types_gene_sizes_intersect$X7
GC_SV_sizes_df <- data.frame(gene, gene_length, GC_group)
#View(GC_SV_sizes_df)
GC_SV_sizes_df <- GC_SV_sizes_df %>% group_by(gene, GC_group) %>% summarise(gene_length=sum(gene_length))
GC_SV_sizes_df$gene_length <- log(GC_SV_sizes_df$gene_length)
library("RColorBrewer")
pdf("ITX_GC_SV_sizes_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~GC_group, data = GC_SV_sizes_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
summary(aov(gene_length~GC_group, data = GC_SV_sizes_df))
table(GC_SV_sizes_df$GC_group)
ITX_GC_stat <- data.frame(GC_SV_sizes_df) %>% wilcox_test(gene_length~GC_group)
write.xlsx(ITX_GC_stat, file="ITX_GC_stat.xlsx")

#TANDEM_DUP
TANDEM_DUP_GC_SV_merged_df <- GC_SV_merged_df %>% filter(SV_type=="TANDEM_DUP")
TANDEM_DUP_GC_SV_intestinal <- TANDEM_DUP_GC_SV_merged_df %>% filter(Group=="intestinal")
write.table(TANDEM_DUP_GC_SV_intestinal, file="TANDEM_DUP_GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
TANDEM_DUP_GC_SV_diffuse <- TANDEM_DUP_GC_SV_merged_df %>% filter(Group=="diffuse")
write.table(TANDEM_DUP_GC_SV_diffuse, file="TANDEM_DUP_GC_SV_diffuse.bed", sep="\t", quote = F, row.names = F, col.names = F)
TANDEM_DUP_GC_SV_mixed <- TANDEM_DUP_GC_SV_merged_df %>% filter(Group=="mixed")
write.table(TANDEM_DUP_GC_SV_mixed, file="TANDEM_DUP_GC_SV_mixed.bed", sep="\t", quote = F, row.names = F, col.names = F)

TANDEM_DUP_GC_merged_types_gene_sizes_intersect <- readr::read_tsv("TANDEM_DUP_GC_merged_types_gene_sizes_intersect.bed", col_names = F) %>% distinct() %>% filter(X8!=0)
#table(GC_merged_types_gene_sizes_intersect$X8)
nrow(TANDEM_DUP_GC_merged_types_gene_sizes_intersect)
gene_length=TANDEM_DUP_GC_merged_types_gene_sizes_intersect$X3-TANDEM_DUP_GC_merged_types_gene_sizes_intersect$X2
gene <- TANDEM_DUP_GC_merged_types_gene_sizes_intersect$X4
GC_group <- TANDEM_DUP_GC_merged_types_gene_sizes_intersect$X7
GC_SV_sizes_df <- data.frame(gene, gene_length, GC_group)
#View(GC_SV_sizes_df)
GC_SV_sizes_df <- GC_SV_sizes_df %>% group_by(gene, GC_group) %>% summarise(gene_length=sum(gene_length))
GC_SV_sizes_df$gene_length <- log(GC_SV_sizes_df$gene_length)
library("RColorBrewer")
pdf("TANDEM_DUP_GC_SV_sizes_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~GC_group, data = GC_SV_sizes_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
summary(aov(gene_length~GC_group, data = GC_SV_sizes_df))
table(GC_SV_sizes_df$GC_group)
TANDEM_DUP_GC_stat <- data.frame(GC_SV_sizes_df) %>% wilcox_test(gene_length~GC_group)
write.xlsx(TANDEM_DUP_GC_stat, file="TANDEM_DUP_GC_stat.xlsx")
