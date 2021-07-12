setwd("C:/Users/user/Downloads/COSMIC")
library(dplyr)
library(openxlsx)

#MUTATIONS
cosmicData <- readr::read_tsv("CosmicCompleteTargetedScreensMutantExport.tsv",col_names = TRUE)
cancer_cosmicData <- cosmicData %>% filter(`Primary histology`==c("carcinoma")) %>% 
  select(Gene=`Gene name`, Gene_length=`Gene CDS length`, Cancer=`Primary site`, Mutation_type=`Mutation Description`, Mutation_CDS=`Mutation CDS`, SNP,
         Coordinate=`Mutation genome position`, Strand=`Mutation strand`) %>% filter(!is.na(Coordinate)) %>% 
  filter(!is.na(Mutation_type), Mutation_type!="Unknown")
#save(Gastric_cancer_cosmicData, file="Gastric_cancer_cosmicData.Rdata")
View(cancer_cosmicData)
chrom <- paste0("chr", sapply(strsplit(as.character(cancer_cosmicData$Coordinate), ":"), "[", 1))
Coordinates <- sapply(strsplit(as.character(cancer_cosmicData$Coordinate), ":"), "[", 2)
chromStart <- as.numeric(sapply(strsplit(as.character(Coordinates), "-"), "[", 1))
chromEnd <- as.numeric(sapply(strsplit(as.character(Coordinates), "-"), "[", 2))
Cancer <- cancer_cosmicData$Cancer 
Mutation <- cancer_cosmicData$Mutation_type
Srtand <- cancer_cosmicData$Strand

#Mutations preparation: universal MUT type
#MUT_strand <- rep("+",length(Srtand))
Mutation[grep("Substitution", Mutation)] <- "Substitution"
Mutation[grep("Deletion", Mutation,ignore.case=TRUE)] <- "Deletion"
Mutation[grep("Insertion", Mutation, ignore.case=TRUE)] <- "Insertion"
table(Mutation)

MUT_df <- data.frame(cbind(chrom, chromStart, chromEnd, Mutation, Cancer, Srtand)) %>% 
  filter(Mutation%in%c("Deletion","Substitution","Insertion"))
dim(MUT_df)
View(MUT_df)
GC_MUT_cosmic_bed <- MUT_df %>% filter(Cancer=="stomach")
write.table(GC_MUT_cosmic_bed, file="GC_MUT_cosmic.bed", sep="\t", quote = F, row.names = F, col.names = F)
OTHER_cancers_MUT_cosmic_bed <- MUT_df %>% filter(Cancer!=c("NS", "stomach"), chromStart<chromEnd)
write.table(OTHER_cancers_MUT_cosmic_bed, file="OTHER_cancers_MUT_cosmic.bed", sep="\t", quote = F, row.names = F, col.names = F)


#Structural variants
SV_GC <- openxlsx::read.xlsx("GastriCancer_Mutations_SV.xlsx", sheet=3)
SV_GC$Chr <- paste0("chr", SV_GC$Chr)
SV_GC_wang2014 <- SV_GC[,c(1:3,5,4,6)]
SV_GC_wang2014$Strand <- rep("+", nrow(SV_GC_wang2014))
write.table(SV_GC_wang2014, file="SV_GC_wang2014.bed", sep="\t", quote = F, row.names = F, col.names = F)

#SV COSMIC
CosmicBreakpointsExport <- readr::read_tsv("CosmicBreakpointsExport.tsv",col_names = TRUE)
# View(CosmicBreakpointsExport)
# unique(CosmicBreakpointsExport$`Primary histology`)
SV_OTHER_cancers_cosmic <- CosmicBreakpointsExport %>% filter(`Primary histology`==c("carcinoma"))
Other_SV_carcinoma_From <- SV_OTHER_cancers_cosmic %>% select(Cancer_type=`Primary site`, Chr=`Chrom From`,Start=`Location From min`,Stop=`Location From max`, SV=`Mutation Type`, Strand=`Strand From`)
Other_SV_carcinoma_To <- SV_OTHER_cancers_cosmic %>% select(Cancer_type=`Primary site`, Chr=`Chrom To`,Start=`Location To min`,Stop=`Location To max`, SV=`Mutation Type`, Strand=`Strand To`)
Other_SV_carcinoma.bed <- rbind(Other_SV_carcinoma_From, Other_SV_carcinoma_To)
unique(Other_SV_carcinoma.bed$SV)
Other_SV_carcinoma.bed$Chr <- paste0("chr",Other_SV_carcinoma.bed$Chr)
#unique(Other_SV_carcinoma.bed$SV)
SV_OTHER_cancers_cosmic_bed <- data.frame(Other_SV_carcinoma.bed[,c(2:4,5,1,6)])
SV_OTHER_cancers_cosmic_bed$Strand <- rep("+",nrow(SV_OTHER_cancers_cosmic_bed))
write.table(SV_OTHER_cancers_cosmic_bed, file="SV_OTHER_cancers_cosmic.bed", sep="\t", quote = F, row.names = F, col.names = F)
#table(SV_OTHER_cancers_cosmic_bed$Strand)

#SV intersect with genes
SV_gene_sizes_intersect <- na.omit(readr::read_tsv("SV_gene_sizes_intersect.bed", col_names = F))
View(SV_gene_sizes_intersect)
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
SV_group <- ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]==0, "GC only", 
                   ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]>0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Shared", 
                          ifelse(SV_gene_sizes_intersect$SV_per_gene[GC_ind]==0 & SV_gene_sizes_intersect$SV_per_gene[Other_ind]>0, "Other only", "Zero")))
table(SV_group)
SV_gene_sizes_grouped_df <- data.frame(SV_gene_sizes_intersect[duplicated(SV_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
SV_gene_sizes_grouped_df$SV_group <- SV_group
SV_gene_sizes_grouped_df <- SV_gene_sizes_grouped_df %>% filter(SV_group!="Zero")
SV_gene_sizes_grouped_df$gene_length <- log(SV_gene_sizes_grouped_df$gene_length)
library("RColorBrewer")
pdf("SV_gene_sizes_grouped_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~SV_group, data = SV_gene_sizes_grouped_df, font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
SV_stat <- wilcox_test(gene_length~SV_group, data = SV_gene_sizes_grouped_df)
write.xlsx(SV_stat, file="SV_stat.Wilcox.test.xlsx")


#MUT intersect with genes
MUT_gene_sizes_intersect <- na.omit(readr::read_tsv("MUT_gene_sizes_intersect.bed", col_names = F)) 
MUT_gene_sizes_intersect$X1[grep("chrX", MUT_gene_sizes_intersect$X1, ignore.case = T)] <- "chrX"
MUT_gene_sizes_intersect <- MUT_gene_sizes_intersect %>% distinct()
View(MUT_gene_sizes_intersect)
#table(MUT_gene_sizes_intersect$X8)
gene <- MUT_gene_sizes_intersect$X4
gene_length <- MUT_gene_sizes_intersect$X3-MUT_gene_sizes_intersect$X2
Cancer_group <- MUT_gene_sizes_intersect$X7
MUT_per_gene <- MUT_gene_sizes_intersect$X8
MUT_gene_sizes_intersect <- data.frame(cbind(gene, gene_length, Cancer_group, MUT_per_gene))
MUT_gene_sizes_intersect$gene_length <- as.numeric(MUT_gene_sizes_intersect$gene_length)
MUT_gene_sizes_intersect$MUT_per_gene <- as.numeric(MUT_gene_sizes_intersect$MUT_per_gene)
str(MUT_gene_sizes_intersect)
#View(MUT_gene_sizes_intersect)
MUT_gene_sizes_intersect <- MUT_gene_sizes_intersect %>% group_by(gene, Cancer_group) %>% summarise(gene_length=sum(gene_length), MUT_per_gene=sum(MUT_per_gene))
View(MUT_gene_sizes_intersect)
GC_ind <- seq(1,nrow(MUT_gene_sizes_intersect),2)
Other_ind <- seq(2,nrow(MUT_gene_sizes_intersect),2)
MUT_group <- ifelse(MUT_gene_sizes_intersect$MUT_per_gene[GC_ind]>0 & MUT_gene_sizes_intersect$MUT_per_gene[Other_ind]==0, "GC only", 
                   ifelse(MUT_gene_sizes_intersect$MUT_per_gene[GC_ind]>0 & MUT_gene_sizes_intersect$MUT_per_gene[Other_ind]>0, "Shared", 
                          ifelse(MUT_gene_sizes_intersect$MUT_per_gene[GC_ind]==0 & MUT_gene_sizes_intersect$MUT_per_gene[Other_ind]>0, "Other only", "Zero")))
table(MUT_group)
MUT_gene_sizes_grouped_df <- data.frame(MUT_gene_sizes_intersect[duplicated(MUT_gene_sizes_intersect$gene),c(1,3)], row.names = NULL)
MUT_gene_sizes_grouped_df$MUT_group <- MUT_group
MUT_gene_sizes_grouped_df <- MUT_gene_sizes_grouped_df %>% filter(MUT_group!="Zero")
MUT_gene_sizes_grouped_df$gene_length <- log(MUT_gene_sizes_grouped_df$gene_length)
library("RColorBrewer")
pdf("MUT_gene_sizes_grouped_df.pdf", width = 5, height = 5.5)
par(font.axis = 2)
boxplot(gene_length~MUT_group, data = MUT_gene_sizes_grouped_df, ylim=c(4,16),font=2, pch=20, cex.axis=1.2, las=1, xlab="", ylab="", col=brewer.pal(n = 3, name = "RdBu"))
dev.off()
library(rstatix)
MUT_stat <- wilcox_test(gene_length~MUT_group, data = MUT_gene_sizes_grouped_df)
write.xlsx(MUT_stat, file="MUT_stat.Wilcox.test.xlsx")


#Check how many mutations/SV per gene?
#SV
SV_per_gene_other_cancers <- na.omit(readr::read_tsv("SV_OTHER_cancers_cosmic.genic", col_names = F))
Cancer <- SV_per_gene_other_cancers$X5
Count_per_gene <- SV_per_gene_other_cancers$X7
Cancer_group <- rep("Other cancers", nrow(SV_per_gene_other_cancers))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
SV_group <- SV_per_gene_other_cancers$X4
SV_freq_other_cancers <- data.frame(Cancer, Cancer_group, SV_group, Count_per_gene, Genic_group)

SV_per_gene_GC <- na.omit(readr::read_tsv("SV_GC_wang2014.genic", col_names = F))
Cancer <- SV_per_gene_GC$X5
Count_per_gene <- SV_per_gene_GC$X7
Cancer_group <- rep("Gastric cancer", nrow(SV_per_gene_GC))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
SV_group <- SV_per_gene_GC$X4
SV_freq_GC <- data.frame(Cancer, Cancer_group, SV_group, Count_per_gene, Genic_group)

SV_per_gene_Random <- na.omit(readr::read_tsv("SV_Random_sites.genic", col_names = F))
Cancer <- SV_per_gene_Random$X5
Count_per_gene <- SV_per_gene_Random$X7
Cancer_group <- rep("Random sites", nrow(SV_per_gene_GC))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
SV_group <- SV_per_gene_Random$X4
SV_freq_Random <- data.frame(Cancer, Cancer_group, SV_group, Count_per_gene, Genic_group)
SV_freq_per_gene <- rbind(SV_freq_GC,SV_freq_other_cancers,SV_freq_Random)
View(SV_freq_per_gene)
save(SV_freq_per_gene, file="SV_freq_per_gene.Rdata")

#MUT
MUT_per_gene_other_cancers <- na.omit(readr::read_tsv("OTHER_cancers_MUT_cosmic.genic", col_names = F))
Cancer <- MUT_per_gene_other_cancers$X5
Count_per_gene <- MUT_per_gene_other_cancers$X7
Cancer_group <- rep("Other cancers", nrow(MUT_per_gene_other_cancers))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
MUT_group <- MUT_per_gene_other_cancers$X4
MUT_freq_other_cancers <- data.frame(Cancer, Cancer_group, MUT_group, Count_per_gene, Genic_group)

MUT_per_gene_GC <- na.omit(readr::read_tsv("GC_MUT_cosmic.genic", col_names = F))
Cancer <- MUT_per_gene_GC$X5
Count_per_gene <- MUT_per_gene_GC$X7
Cancer_group <- rep("Gastric cancer", nrow(MUT_per_gene_GC))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
MUT_group <- MUT_per_gene_GC$X4
MUT_freq_GC <- data.frame(Cancer, Cancer_group, MUT_group, Count_per_gene, Genic_group)

MUT_per_gene_Random <- na.omit(readr::read_tsv("MUT_Random_sites.genic", col_names = F))
Cancer <- MUT_per_gene_Random$X5
Count_per_gene <- MUT_per_gene_Random$X7
Cancer_group <- rep("Random sites", nrow(MUT_per_gene_GC))
Genic_group <- ifelse(Count_per_gene>0, "genic", "Intergenic")
MUT_group <- MUT_per_gene_Random$X4
MUT_freq_Random <- data.frame(Cancer, Cancer_group, MUT_group, Count_per_gene, Genic_group)
MUT_freq_per_gene <- rbind(MUT_freq_GC,MUT_freq_other_cancers,MUT_freq_Random)
View(MUT_freq_per_gene)
save(MUT_freq_per_gene, file="MUT_freq_per_gene.Rdata")
