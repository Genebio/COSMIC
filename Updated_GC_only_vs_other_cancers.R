setwd("C:/Users/user/Downloads/COSMIC")
#Structural variants
SV_GC <- openxlsx::read.xlsx("GastriCancer_Mutations_SV.xlsx", sheet=3)
SV_GC$Chr <- paste0("chr", SV_GC$Chr)
SV_GC_wang2014 <- SV_GC[,c(1:3,5,4,6)]
SV_GC_wang2014$Strand <- rep("+", nrow(SV_GC_wang2014))
write.table(SV_GC_wang2014, file="SV_GC_wang2014.bed", sep="\t", quote = F, row.names = F, col.names = F)

#SV GC
load("GC_SV_merged_df.Rdata")
str(GC_SV_merged_df)
names(GC_SV_merged_df) <- c("Chr", "Start", "Stop", "SV", "Cancer_type", "Strand")
GC_SV_merged_df$Cancer_type <- rep("GC", nrow(GC_SV_merged_df))
GC_SV <- GC_SV_merged_df
#write.table(GC_SV_bed, file="GC_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
merged_all <- rbind(GC_SV_merged_df,SV_OTHER_cancers_cosmic_bed) %>% distinct()
View(table(merged_all$Cancer_type))

#SV
breast_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="breast")
write.table(breast_SV, file="breast_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
ovary_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="ovary")
write.table(ovary_SV, file="ovary_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
# pancreas_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="pancreas")
# write.table(pancreas_SV, file="pancreas_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
prostate_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="prostate")
write.table(prostate_SV, file="prostate_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)

#GC vs other
GC_other <- rbind(breast_SV,ovary_SV,prostate_SV)
write.table(GC_other, file="GC_other.bed", sep="\t", quote = F, row.names = F, col.names = F)
breast_other <- rbind(GC_SV,ovary_SV,prostate_SV)
write.table(breast_other, file="breast_other.bed", sep="\t", quote = F, row.names = F, col.names = F)
ovary_other <- rbind(GC_SV,breast_SV,prostate_SV)
write.table(ovary_other, file="ovary_other.bed", sep="\t", quote = F, row.names = F, col.names = F)
prostate_other <- rbind(GC_SV,breast_SV,ovary_SV)
write.table(prostate_other, file="prostate_other.bed", sep="\t", quote = F, row.names = F, col.names = F)




#GC intersect with genes
SV_gene_sizes_intersect <- na.omit(readr::read_tsv("GC_other_gene_sizes_intersect.bed", col_names = F))
#View(SV_gene_sizes_intersect)
gene <- SV_gene_sizes_intersect$X4
gene_length <- SV_gene_sizes_intersect$X3-SV_gene_sizes_intersect$X2
Cancer_group <- SV_gene_sizes_intersect$X7
SV_per_gene <- SV_gene_sizes_intersect$X8
SV_gene_sizes_intersect <- data.frame(cbind(gene, gene_length, Cancer_group, SV_per_gene))
SV_gene_sizes_intersect$gene_length <- as.numeric(SV_gene_sizes_intersect$gene_length)
SV_gene_sizes_intersect$SV_per_gene <- as.numeric(SV_gene_sizes_intersect$SV_per_gene)
str(SV_gene_sizes_intersect)
#View(SV_gene_sizes_intersect)
SV_gene_sizes_intersect <- SV_gene_sizes_intersect %>% group_by(gene, Cancer_group) %>% summarise(gene_length=sum(gene_length), SV_per_gene=sum(SV_per_gene))
head(SV_gene_sizes_intersect)
GC_ind <- seq(1,nrow(SV_gene_sizes_intersect),2)
Other_ind <- seq(2,nrow(SV_gene_sizes_intersect),2)
SV_group <- ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]==0, "GC only", 
                   ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Shared", 
                          ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]==0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Other only", "Zero")))
table(SV_group)
SV_gene_sizes_grouped_df <- data.frame(SV_gene_sizes_intersect[duplicated(SV_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
SV_gene_sizes_grouped_df$SV_group <- SV_group
SV_gene_sizes_grouped_df <- SV_gene_sizes_grouped_df %>% filter(SV_group!="Zero")
SV_gene_sizes_grouped_df$gene_length <- log(SV_gene_sizes_grouped_df$gene_length)
GC_sizes_genes <- SV_gene_sizes_grouped_df %>% filter(SV_group=="GC only")
library("RColorBrewer")
pdf("GC_other.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~SV_group, data = SV_gene_sizes_grouped_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
GC_other_stat_df <- SV_gene_sizes_grouped_df
save(GC_other_stat_df, file="GC_other_stat_df.Rdata")
GC_other_stat <- t_test(gene_length~SV_group, data = SV_gene_sizes_grouped_df)
write.xlsx(GC_other_stat, file="GC_other_stat.xlsx")


#breast intersect with genes
SV_gene_sizes_intersect <- na.omit(readr::read_tsv("breast_other_gene_sizes_intersect.bed", col_names = F))
#View(SV_gene_sizes_intersect)
gene <- SV_gene_sizes_intersect$X4
gene_length <- SV_gene_sizes_intersect$X3-SV_gene_sizes_intersect$X2
Cancer_group <- SV_gene_sizes_intersect$X7
SV_per_gene <- SV_gene_sizes_intersect$X8
SV_gene_sizes_intersect <- data.frame(cbind(gene, gene_length, Cancer_group, SV_per_gene))
SV_gene_sizes_intersect$gene_length <- as.numeric(SV_gene_sizes_intersect$gene_length)
SV_gene_sizes_intersect$SV_per_gene <- as.numeric(SV_gene_sizes_intersect$SV_per_gene)
str(SV_gene_sizes_intersect)
#View(SV_gene_sizes_intersect)
SV_gene_sizes_intersect <- SV_gene_sizes_intersect %>% group_by(gene, Cancer_group) %>% summarise(gene_length=sum(gene_length), SV_per_gene=sum(SV_per_gene))
GC_ind <- seq(1,nrow(SV_gene_sizes_intersect),2)
Other_ind <- seq(2,nrow(SV_gene_sizes_intersect),2)
SV_group <- ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]==0, "Breast only", 
                   ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Shared", 
                          ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]==0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Other only", "Zero")))
SV_gene_sizes_grouped_df <- data.frame(SV_gene_sizes_intersect[duplicated(SV_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
SV_gene_sizes_grouped_df$SV_group <- SV_group
SV_gene_sizes_grouped_df <- SV_gene_sizes_grouped_df %>% filter(SV_group!="Zero")
SV_gene_sizes_grouped_df$gene_length <- log(SV_gene_sizes_grouped_df$gene_length)
breast_sizes_genes <- SV_gene_sizes_grouped_df %>% filter(SV_group=="Breast only")
library("RColorBrewer")
pdf("breast_other.pdf", width = 5.1, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~SV_group, data = SV_gene_sizes_grouped_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
breast_other_stat_df <- SV_gene_sizes_grouped_df
save(breast_other_stat_df, file="breast_other_stat_df.Rdata")
breast_other_stat <- wilcox_test(gene_length~SV_group, data = SV_gene_sizes_grouped_df)
write.xlsx(breast_other_stat, file="breast_other_stat.xlsx")
breast_other_stat

#ovary intersect with genes
SV_gene_sizes_intersect <- na.omit(readr::read_tsv("ovary_other_gene_sizes_intersect.bed", col_names = F))
#View(SV_gene_sizes_intersect)
gene <- SV_gene_sizes_intersect$X4
gene_length <- SV_gene_sizes_intersect$X3-SV_gene_sizes_intersect$X2
Cancer_group <- SV_gene_sizes_intersect$X7
SV_per_gene <- SV_gene_sizes_intersect$X8
SV_gene_sizes_intersect <- data.frame(cbind(gene, gene_length, Cancer_group, SV_per_gene))
SV_gene_sizes_intersect$gene_length <- as.numeric(SV_gene_sizes_intersect$gene_length)
SV_gene_sizes_intersect$SV_per_gene <- as.numeric(SV_gene_sizes_intersect$SV_per_gene)
str(SV_gene_sizes_intersect)
#View(SV_gene_sizes_intersect)
SV_gene_sizes_intersect <- SV_gene_sizes_intersect %>% group_by(gene, Cancer_group) %>% summarise(gene_length=sum(gene_length), SV_per_gene=sum(SV_per_gene))
GC_ind <- seq(2,nrow(SV_gene_sizes_intersect),2)
Other_ind <- seq(1,nrow(SV_gene_sizes_intersect),2)
SV_group <- ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]==0, "Ovary only", 
                   ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Shared", 
                          ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]==0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Other only", "Zero")))
#SV_group <- factor(SV_group, levels = c("Ovary only","Other only","Shared","Zero"))
table(SV_group)
SV_gene_sizes_grouped_df <- data.frame(SV_gene_sizes_intersect[duplicated(SV_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
SV_gene_sizes_grouped_df$SV_group <- SV_group
SV_gene_sizes_grouped_df <- SV_gene_sizes_grouped_df %>% filter(SV_group!="Zero")
SV_gene_sizes_grouped_df$gene_length <- log(SV_gene_sizes_grouped_df$gene_length)
SV_gene_sizes_grouped_df$SV_group <- factor(SV_gene_sizes_grouped_df$SV_group, levels = c("Ovary only","Other only","Shared"))
ovary_sizes_genes <- SV_gene_sizes_grouped_df %>% filter(SV_group=="Ovary only")
library("RColorBrewer")
pdf("ovary_other.pdf", width = 5.1, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~SV_group, data = SV_gene_sizes_grouped_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
ovary_other_stat_df <- SV_gene_sizes_grouped_df
save(ovary_other_stat_df, file="ovary_other_stat_df.Rdata")
ovary_other_stat <- wilcox_test(gene_length~SV_group, data = SV_gene_sizes_grouped_df)
write.xlsx(ovary_other_stat, file="ovary_other_stat.xlsx")

#prostate intersect with genes
SV_gene_sizes_intersect <- na.omit(readr::read_tsv("prostate_other_gene_sizes_intersect.bed", col_names = F))
#View(SV_gene_sizes_intersect)
gene <- SV_gene_sizes_intersect$X4
gene_length <- SV_gene_sizes_intersect$X3-SV_gene_sizes_intersect$X2
Cancer_group <- SV_gene_sizes_intersect$X7
SV_per_gene <- SV_gene_sizes_intersect$X8
SV_gene_sizes_intersect <- data.frame(cbind(gene, gene_length, Cancer_group, SV_per_gene))
SV_gene_sizes_intersect$gene_length <- as.numeric(SV_gene_sizes_intersect$gene_length)
SV_gene_sizes_intersect$SV_per_gene <- as.numeric(SV_gene_sizes_intersect$SV_per_gene)
str(SV_gene_sizes_intersect)
#View(SV_gene_sizes_intersect)
SV_gene_sizes_intersect <- SV_gene_sizes_intersect %>% group_by(gene, Cancer_group) %>% summarise(gene_length=sum(gene_length), SV_per_gene=sum(SV_per_gene))
head(SV_gene_sizes_intersect)
GC_ind <- seq(2,nrow(SV_gene_sizes_intersect),2)
Other_ind <- seq(1,nrow(SV_gene_sizes_intersect),2)
SV_group <- ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]==0, "Prostate only", 
                   ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Shared", 
                          ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]==0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Other only", "Zero")))
#SV_group <- factor(SV_group, levels = c("Prostate only","Other only","Shared","Zero"))
table(SV_group)
SV_gene_sizes_grouped_df <- data.frame(SV_gene_sizes_intersect[duplicated(SV_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
SV_gene_sizes_grouped_df$SV_group <- SV_group
SV_gene_sizes_grouped_df <- SV_gene_sizes_grouped_df %>% filter(SV_group!="Zero")
SV_gene_sizes_grouped_df$gene_length <- log(SV_gene_sizes_grouped_df$gene_length)
SV_gene_sizes_grouped_df$SV_group <- factor(SV_gene_sizes_grouped_df$SV_group, levels = c("Prostate only","Other only","Shared"))
prostate_sizes_genes <- SV_gene_sizes_grouped_df %>% filter(SV_group=="Prostate only")
library("RColorBrewer")
pdf("prostate_other.pdf", width = 5.5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~SV_group, data = SV_gene_sizes_grouped_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
prostate_other_stat_df <- SV_gene_sizes_grouped_df
save(prostate_other_stat_df, file="prostate_other_stat_df.Rdata")
prostate_other_stat <- wilcox_test(gene_length~SV_group, data = SV_gene_sizes_grouped_df)
write.xlsx(prostate_other_stat, file="prostate_other_stat.xlsx")

#Data frame of 'Only" genes - not shared
Merged_genes_sizes <- rbind(GC_sizes_genes, breast_sizes_genes, ovary_sizes_genes, prostate_sizes_genes)
save(Merged_genes_sizes, file="Merged_genes_sizes.Rdata")

#Merge the data with transcriprion data
RNAseq_data <- openxlsx::read.xlsx("AGS_RANAseq_transcription.xlsx")
head(RNAseq_data)
gene_length <- RNAseq_data$end-RNAseq_data$start
RNAseq_genes_transcription <- data.frame(gene=RNAseq_data$gene_name, gene_length=log(gene_length), transcription=log(rowMeans(RNAseq_data[,9:10]+1)))
RNAseq_genes_transcription <- RNAseq_genes_transcription[na.omit(match(toupper(SV_gene_sizes_grouped_df$gene), toupper(RNAseq_genes_transcription$gene))),]
Group <- SV_gene_sizes_grouped_df$SV_group[SV_gene_sizes_grouped_df$gene%in%RNAseq_data$gene_name]
RNAseq_genes_transcription$Group <- Group
table(RNAseq_genes_transcription$Group)

# RNAseq_data <- openxlsx::read.xlsx("AGS_RANAseq_transcription.xlsx")
# head(RNAseq_data)
# RNAseq_genes_transcription <- data.frame(gene=RNAseq_data$gene_name, transcription=log(rowMeans(RNAseq_data[,9:10]+1)))
# na.omit(match(SV_gene_sizes_grouped_df$gene, RNAseq_genes_transcription$gene))
# Transcription <- RNAseq_genes_transcription$transcription[na.omit(match(toupper(SV_gene_sizes_grouped_df$gene), toupper(RNAseq_genes_transcription$gene)))]
# Merged_genes_sizes_matched <- SV_gene_sizes_grouped_df[SV_gene_sizes_grouped_df$gene%in%RNAseq_data$gene_name,]
# Merged_genes_sizes_matched$Transcription <- Transcription
# Merged_genes_sizes_transcription <- Merged_genes_sizes_matched
# save(Merged_genes_sizes_transcription, file="Merged_genes_sizes_transcription.Rdata")
# table(Merged_genes_sizes_transcription$SV_group)

# Merged_genes_sizes_transcription$SV_group <- ifelse(Merged_genes_sizes_transcription$SV_group=="GC only", "GC only", "Other")
# #Make a nice scatterplot of gene_length to expression by cancer type
# SV_gene_sizes_grouped_df$SV_group <- ifelse(SV_gene_sizes_grouped_df$SV_group=="GC only", "GC only", "Other")
x=RNAseq_genes_transcription
head(x)
gcind=which(x[,4]=="GC only")
ind=c(gcind,sample(which(x[,4]=="Other"),length(gcind)))
data=x[ind,c(2,3)]
label=x[ind,4]
cols=rep("black",length(label))
cols[label=="GC only"]="red3"
pdf("GC_only_other_transcription.pdf")
plot(data[,1],data[,2],col=cols,pch=20,xlab="Gene Size (log(bp))",ylab="Transcription (log(#reads))",cex.lab=1.4,cex.axis=1.4)
dev.off()