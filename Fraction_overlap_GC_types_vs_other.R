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

#Check chromosomes in REP files

#HG19 reference genes
HG19_genes_check <- readr::read_tsv("HG19.GENES.bed", col_names = F)
unique(HG19_genes_check$X1)
for (i in 1:23){
  if (any(grepl(paste0("chr",i,"_"),HG19_genes_check$X1))){
    HG19_genes_check$X1[grep(paste0("chr",i,"_"),HG19_genes_check$X1)] <- paste0("chr",i)
  }}
HG19_genes_check$X1[grep("chrUn",HG19_genes_check$X1)] <- "chrM"
HG19_genes_check$X1[grep("chrX",HG19_genes_check$X1)] <- "chrX"
HG19_genes_check$X1[grep("chrY",HG19_genes_check$X1)] <- "chrY"
HG19_genes_check <- HG19_genes_check %>% filter(!X1%in%c("chrX","chrY","chrM","chrMT"))
unique(HG19_genes_check$X1)
write.table(HG19_genes_check, file="HG19_genes.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)

#REP
REP_check <- readr::read_tsv("REP.HG19.bed", col_names = F)
unique(REP_check$X1)
for (i in 1:23){
  if (any(grepl(paste0("chr",i,"_"),REP_check$X1))){
    REP_check$X1[grep(paste0("chr",i,"_"),REP_check$X1)] <- paste0("chr",i)
  }}
REP_check$X1[grep("chrUn",REP_check$X1)] <- "chrM"
REP_check$X1[grep("chrX",REP_check$X1)] <- "chrX"
REP_check$X1[grep("chrY",REP_check$X1)] <- "chrY"
REP_check <- REP_check %>% filter(!X1%in%c("chrX","chrY","chrM","chrMT"))
unique(REP_check$X1)
write.table(REP_check, file="REP.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)

#polyA
polyA_check <- readr::read_tsv("Arep.HG19.bed", col_names = F)
unique(polyA_check$X1)
polyA_check <- polyA_check %>% filter(!X1%in%c("chrX","chrY"))
unique(polyA_check$X1)
write.table(polyA_check, file="polyA.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)

#polyTA
polyTA_check <- readr::read_tsv("TArep.HG19.bed", col_names = F)
unique(polyTA_check$X1)
polyTA_check <- polyTA_check %>% filter(!X1%in%c("chrX","chrY"))
unique(polyTA_check$X1)
write.table(polyTA_check, file="polyTA.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)

#TSS
TSS_check <- readr::read_tsv("TSS.HG19.bed", col_names = F)
for (i in 1:23){
  if (any(grepl(paste0("chr",i,"_"),TSS_check$X1))){
    TSS_check$X1[grep(paste0("chr",i,"_"),TSS_check$X1)] <- paste0("chr",i)
  }}
TSS_check$X1[grep("chrUn",TSS_check$X1)] <- "chrM"
TSS_check$X1[grep("chrX",TSS_check$X1)] <- "chrX"
TSS_check$X1[grep("chrY",TSS_check$X1)] <- "chrY"
TSS_check <- TSS_check %>% filter(!X1%in%c("chrX","chrY","chrM"))
unique(TSS_check$X1)
write.table(TSS_check, file="TSS.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)

#TES
TES_check <- readr::read_tsv("TES.HG19.bed", col_names = F)
for (i in 1:23){
  if (any(grepl(paste0("chr",i,"_"),TES_check$X1))){
    TES_check$X1[grep(paste0("chr",i,"_"),TES_check$X1)] <- paste0("chr",i)
  }}
TES_check$X1[grep("chrUn",TES_check$X1)] <- "chrM"
TES_check$X1[grep("chrX",TES_check$X1)] <- "chrX"
TES_check$X1[grep("chrY",TES_check$X1)] <- "chrY"
TES_check <- TES_check %>% filter(!X1%in%c("chrX","chrY","chrM"))
unique(TES_check$X1)
write.table(TES_check, file="TES.HG19_clean.bed", sep="\t", quote = F, row.names = F, col.names = F)
