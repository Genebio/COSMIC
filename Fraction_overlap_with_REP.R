#Find fraction overlap with repetitive elements
REP <- readr::read_tsv("REP_FLANK.5000_Int_other_SV.bed", col_names = F)
View(REP)
gcind <- seq(1,nrow(REP),2)
other_ind <- seq(2,nrow(REP),2)
gc_overlaped <- round(sum(REP$X8[gcind]>0)*200/length(REP$X8))
gc_overlaped
other_overlaped <- round(sum(REP$X8[other_ind]>0)*200/length(REP$X8))
other_overlaped
group <- ifelse(REP$X8[gcind]>0 & REP$X8[other_ind]==0, "GC only", 
                   ifelse(REP$X8[gcind]>0 & REP$X8[other_ind]>0, "Shared", 
                          ifelse(REP$X8[gcind]==0 & REP$X8[other_ind]>0, "Other only", "Zero")))
#group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

Arep.HG19 <- readr::read_tsv("Arep.HG19_FLANK.5000_Int_other_SV.bed", col_names = F)
gcind <- seq(1,nrow(Arep.HG19),2)
other_ind <- seq(2,nrow(Arep.HG19),2)
gc_overlaped <- round(sum(Arep.HG19$X8[gcind]>0)*200/length(Arep.HG19$X8))
gc_overlaped
other_overlaped <- round(sum(Arep.HG19$X8[other_ind]>0)*200/length(Arep.HG19$X8))
other_overlaped
group <- ifelse(Arep.HG19$X8[gcind]>0 & Arep.HG19$X8[other_ind]==0, "GC only", 
                ifelse(Arep.HG19$X8[gcind]>0 & Arep.HG19$X8[other_ind]>0, "Shared", 
                       ifelse(Arep.HG19$X8[gcind]==0 & Arep.HG19$X8[other_ind]>0, "Other only", "Zero")))
#group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

TArep.HG19 <- readr::read_tsv("TArep.HG19_FLANK.5000_Int_other_SV.bed", col_names = F)
gcind <- seq(1,nrow(TArep.HG19),2)
other_ind <- seq(2,nrow(TArep.HG19),2)
gc_overlaped <- round(sum(TArep.HG19$X8[gcind]>0)*200/length(TArep.HG19$X8))
gc_overlaped
other_overlaped <- round(sum(TArep.HG19$X8[other_ind]>0)*200/length(TArep.HG19$X8))
other_overlaped
group <- ifelse(TArep.HG19$X8[gcind]>0 & TArep.HG19$X8[other_ind]==0, "GC only", 
                ifelse(TArep.HG19$X8[gcind]>0 & TArep.HG19$X8[other_ind]>0, "Shared", 
                       ifelse(TArep.HG19$X8[gcind]==0 & TArep.HG19$X8[other_ind]>0, "Other only", "Zero")))
#group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

TSS.HG19 <- readr::read_tsv("TSS.HG19_FLANK.5000_Int_other_SV.bed", col_names = F)
gcind <- seq(1,nrow(TSS.HG19),2)
other_ind <- seq(2,nrow(TSS.HG19),2)
gc_overlaped <- round(sum(TSS.HG19$X8[gcind]>0)*200/length(TSS.HG19$X8))
gc_overlaped
other_overlaped <- round(sum(TSS.HG19$X8[other_ind]>0)*200/length(TSS.HG19$X8))
other_overlaped
group <- ifelse(TSS.HG19$X8[gcind]>0 & TSS.HG19$X8[other_ind]==0, "GC only", 
                ifelse(TSS.HG19$X8[gcind]>0 & TSS.HG19$X8[other_ind]>0, "Shared", 
                       ifelse(TSS.HG19$X8[gcind]==0 & TSS.HG19$X8[other_ind]>0, "Other only", "Zero")))
#group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

TES.HG19 <- readr::read_tsv("TES.HG19_FLANK.5000_Int_other_SV.bed", col_names = F)
gcind <- seq(1,nrow(TES.HG19),2)
other_ind <- seq(2,nrow(TES.HG19),2)
gc_overlaped <- round(sum(TES.HG19$X8[gcind]>0)*200/length(TES.HG19$X8))
gc_overlaped
other_overlaped <- round(sum(TES.HG19$X8[other_ind]>0)*200/length(TES.HG19$X8))
other_overlaped
group <- ifelse(TES.HG19$X8[gcind]>0 & TES.HG19$X8[other_ind]==0, "GC only", 
                ifelse(TES.HG19$X8[gcind]>0 & TES.HG19$X8[other_ind]>0, "Shared", 
                       ifelse(TES.HG19$X8[gcind]==0 & TES.HG19$X8[other_ind]>0, "Other only", "Zero")))
group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

#Intergenic part
REP <- readr::read_tsv("REP.HG19_FLANK.5000_Int_other_SV_intergenic.bed", col_names = F)
head(REP)
gcind <- seq(1,nrow(REP),2)
other_ind <- seq(2,nrow(REP),2)
gc_overlaped <- round(sum(REP$X8[gcind]>0)*200/length(REP$X8))
gc_overlaped
other_overlaped <- round(sum(REP$X8[other_ind]>0)*200/length(REP$X8))
other_overlaped
group <- ifelse(REP$X8[gcind]>0 & REP$X8[other_ind]==0, "GC only", 
                ifelse(REP$X8[gcind]>0 & REP$X8[other_ind]>0, "Shared", 
                       ifelse(REP$X8[gcind]==0 & REP$X8[other_ind]>0, "Other only", "Zero")))
group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq

Arep <- readr::read_tsv("Arep.HG19_FLANK.5000_Int_other_SV_intergenic.bed", col_names = F)
head(Arep)
gcind <- seq(1,nrow(Arep),2)
other_ind <- seq(2,nrow(Arep),2)
gc_overlaped <- round(sum(Arep$X8[gcind]>0)*200/length(Arep$X8))
gc_overlaped
other_overlaped <- round(sum(Arep$X8[other_ind]>0)*200/length(Arep$X8))
other_overlaped
group <- ifelse(Arep$X8[gcind]>0 & Arep$X8[other_ind]==0, "GC only", 
                ifelse(Arep$X8[gcind]>0 & Arep$X8[other_ind]>0, "Shared", 
                       ifelse(Arep$X8[gcind]==0 & Arep$X8[other_ind]>0, "Other only", "Zero")))
group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)),1)
freq

TArep <- readr::read_tsv("TArep.HG19_FLANK.5000_Int_other_SV_intergenic.bed", col_names = F)
head(TArep)
gcind <- seq(1,nrow(TArep),2)
other_ind <- seq(2,nrow(TArep),2)
gc_overlaped <- round(sum(TArep$X8[gcind]>0)*200/length(TArep$X8))
gc_overlaped
other_overlaped <- round(sum(TArep$X8[other_ind]>0)*200/length(TArep$X8))
other_overlaped
group <- ifelse(TArep$X8[gcind]>0 & TArep$X8[other_ind]==0, "GC only", 
                ifelse(TArep$X8[gcind]>0 & TArep$X8[other_ind]>0, "Shared", 
                       ifelse(TArep$X8[gcind]==0 & TArep$X8[other_ind]>0, "Other only", "Zero")))
group <- group[group!="Zero"]
freq <- round(table(group)*100/sum(table(group)))
freq
