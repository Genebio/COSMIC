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
