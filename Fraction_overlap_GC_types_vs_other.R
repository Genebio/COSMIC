#Fraction overlap with repetitive elements of GC_intestinal, GC_other, prostate, breast, ovary cancers SV
setwd("C:/Users/user/Downloads/COSMIC")
library(dplyr)
library(openxlsx)


clinical_data <- read.xlsx("Updated_GC_SV_wang2014.xlsx", sheet=1)
str(clinical_data)
SV_data <- read.xlsx("Updated_GC_SV_wang2014.xlsx", sheet=2)
str(SV_data)
clinical_data_matched <- data.frame(clinical_data[match(SV_data$Sample, clinical_data$Sample),2:3], row.names=NULL)
chrom <- paste0("chr", SV_data$Chr)
GC_SV_merged_df <- data.frame(chrom, Start=SV_data$Start, End=SV_data$End, Group=clinical_data_matched$GC_type, SV_type=SV_data$SV, Strand=SV_data$Strand)
write.table(GC_SV_merged_df, file="GC_SV_merged_df.bed", sep="\t", quote = F, row.names = F, col.names = F)

#GC_SV_intestinal.bed
GC_SV_intestinal <- GC_SV_merged_df %>% filter(Group=="intestinal",chrom!="chrY")
write.table(GC_SV_intestinal, file="GC_SV_intestinal.bed", sep="\t", quote = F, row.names = F, col.names = F)
unique(GC_SV_intestinal$chrom)
#GC_SV_otherGC.bed
GC_SV_otherGC <- GC_SV_merged_df %>% filter(Group!="intestinal",chrom!="chrY")
write.table(GC_SV_otherGC, file="GC_SV_otherGC.bed", sep="\t", quote = F, row.names = F, col.names = F)
unique(GC_SV_otherGC$chrom)

load("SV_OTHER_cancers_cosmic_bed.Rdata")
#breast_SV.bed
breast_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="breast", Chr!="chr24")
write.table(breast_SV, file="breast_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
unique(breast_SV$Chr)
#ovary_SV.bed
ovary_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="ovary", Chr!="chr24")
write.table(ovary_SV, file="ovary_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
unique(ovary_SV$Chr)
#prostate_SV.bed
prostate_SV <- SV_OTHER_cancers_cosmic_bed %>% filter(Cancer_type=="prostate", Chr!="chr24",Chr!="chr25")
write.table(prostate_SV, file="prostate_SV.bed", sep="\t", quote = F, row.names = F, col.names = F)
unique(prostate_SV$Chr)


