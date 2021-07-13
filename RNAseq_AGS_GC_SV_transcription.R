setwd("C:/Users/user/Downloads/COSMIC")
library(dplyr)
library(openxlsx)
load("SV_gene_sizes_grouped_df.Rdata")
str(SV_gene_sizes_grouped_df)
GC_other_SV_df <- SV_gene_sizes_grouped_df
GC_other_SV_df$SV_group <- ifelse(GC_other_SV_df$SV_group=="GC only", "GC only", "Other")
nrow(GC_other_SV_df)

#Gene transcription level
AGS_RANAseq_transcription <- read.xlsx("AGS_RANAseq_transcription.xlsx")
names(AGS_RANAseq_transcription)
AGS_RANAseq_transcription <- AGS_RANAseq_transcription[na.omit(match(GC_other_SV_df$gene, AGS_RANAseq_transcription$gene_name)),c(7,9:10)]
GC_other_SV_df <- GC_other_SV_df[GC_other_SV_df$gene%in%AGS_RANAseq_transcription$gene_name,]
GC_other_SV_df$log_transcription <- log(rowMeans(AGS_RANAseq_transcription[,2:3]+1))
save(GC_other_SV_df, file="GC_other_SV_df.Rdata")

#Scatterplot of transcription / gene size [1:gene;2:gene_length;3:SV_group;4:log_transcr]
load("GC_other_SV_df.Rdata")
names(x)
x=GC_other_SV_df
gcind=which(x[,3]=="GC only")
ind=c(gcind,sample(which(x[,3]=="Other"),length(gcind)))
data=x[ind,c(2,4)]
label=x[ind,3]
cols=rep("black",length(label))
cols[label=="GC only"]="red3"
pdf("RNAseq_AGS_SV_transcription.pdf")
plot(data[,1],data[,2],col=cols,pch=20,xlab="Gene Size (log(bp))",ylab="Transcription (log(#reads))",cex.lab=1.4,cex.axis=1.4, las=1)
dev.off()