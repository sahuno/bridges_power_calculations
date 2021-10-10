##date Jan 18th 2021
##author; sam Ahuno
###purpose; get biallalleic inactivaction for all tumor types


library(data.table)
library(tidyverse)
library(biomaRt)
library(GenomeInfoDb)
library(GenomicFeatures)
library(reshape2)
library(plyr)
library(ggpubr)
library(RColorBrewer)
library(hablar)
library(ggplot2)



##read the tcga absolute master calls 
df.tcga.mastercalls.abs <- fread("/Users/samuelahuno/Downloads/TCGA_mastercalls.abs_segtabs.fixed.txt", sep = "\t", header= TRUE)
##sanity checks
head(df.tcga.mastercalls.abs)
dim(df.tcga.mastercalls.abs)
###how many unique sample
df.tcga.mastercalls.abs %>% distinct(Sample) %>% dplyr::summarise(n()) #1 11084

###retrive cancer types
df.tcga.aneuploidy <- read.delim("/Users/samuelahuno/Downloads/PANCAN_ArmCallsAndAneuploidyScore_092817.txt", sep = "\t", header= TRUE)
brca.ids <- df.tcga.aneuploidy[which(df.tcga.aneuploidy$Type == "BRCA"),]
##merged table to append tumor types key
merged.tcga.mastercalls.abs.tumor.ID <- merge(df.tcga.mastercalls.abs, df.tcga.aneuploidy,by.x="Sample",by.y="Sample",no.dups=TRUE)
dim(merged.tcga.mastercalls.abs.tumor.ID)
 table(merged.tcga.mastercalls.abs.tumor.ID$Type)
 
 
 #test
 #merge(any.tcga.file.with.sampleIDs, brca.ids,by.x="Sample.x",by.y="Sample",no.dups=TRUE)
 
##get loh samples
#loh.merged.tcga.mastercalls.abs.tumor.ID <- merged.tcga.mastercalls.abs.tumor.ID %>% filter(LOH == 1)


### read breast cacner genes from bridges study
source("/Users/samuelahuno/Polak_lab10082019/bridges_breast/scripts/getCytobandsHg19_bridges_genes.R")
cordinates.bridges.genes.hg19

##find overlaps
gr.loh.merged.tcga.mastercalls.abs.tumor.ID <- makeGRangesFromDataFrame(merged.tcga.mastercalls.abs.tumor.ID,start.field=  "Start", end.field= "End", seqnames.field="Chromosome", keep.extra.columns = T)
ov.gr.loh.merged.tcga.mastercalls.abs.tumor.ID<- findOverlaps(query=cordinates.bridges.genes.hg19, subject=gr.loh.merged.tcga.mastercalls.abs.tumor.ID, type = "within") ##complete overlaps

####append, gene names to columns
gr.loh.merged.tcga.mastercalls.abs.tumor.ID$bridges_genes <- NA
mcols(gr.loh.merged.tcga.mastercalls.abs.tumor.ID)$bridges_genes[subjectHits(ov.gr.loh.merged.tcga.mastercalls.abs.tumor.ID)] <- mcols(cordinates.bridges.genes.hg19)$hgnc_symbol[queryHits(ov.gr.loh.merged.tcga.mastercalls.abs.tumor.ID)]


##convert to data frame
df.loh.merged.tcga.mastercalls.abs.tumor.ID <- as.data.frame(gr.loh.merged.tcga.mastercalls.abs.tumor.ID)
df.loh.merged.tcga.mastercalls.abs.tumor.ID <- df.loh.merged.tcga.mastercalls.abs.tumor.ID %>% drop_na(bridges_genes)
dim(df.loh.merged.tcga.mastercalls.abs.tumor.ID)


###how many patiets are there in each cancer type
total_number_cases_tcga <- df.loh.merged.tcga.mastercalls.abs.tumor.ID %>% distinct_at(vars(Sample),.keep_all=T) %>% group_by(Type) %>% dplyr::summarise(n = n()) ### here breast caner are 28 patients less


##extrcat only loh
stats.loh.tcga.all.tumor.types <- df.loh.merged.tcga.mastercalls.abs.tumor.ID %>% filter(Homozygous_deletion == 1) %>% group_by(Type,bridges_genes) %>% dplyr::summarise(n = n()) #%>% mutate(freq = n /1048,freq_AF = (n /1048)*0.5)
view(stats.loh.tcga.all.tumor.types)
stats.loh.tcga.all.tumor.types %>%filter(Type=="BRCA" & bridges_genes == "BARD1") #%>%

#which patient has homozygous deletions
df.loh.merged.tcga.mastercalls.abs.tumor.ID %>% filter(Homozygous_deletion == 1) %>%filter(Type=="BRCA" & bridges_genes == "BARD1")

stats.loh.tcga.all.tumor.types <- merge(stats.loh.tcga.all.tumor.types,total_number_cases_tcga,by="Type")
setnames(stats.loh.tcga.all.tumor.types,old=c("n.x","n.y"),new = c("ncases_loh",  "total_sample"))




stats.loh.tcga.all.tumor.types.AF <- stats.loh.tcga.all.tumor.types %>% mutate(freq = ncases_loh/total_sample,freq_AF = (ncases_loh/total_sample)*0.5,EI=1/(1-freq_AF))
##prof check 1/(1-0.333333333)
# 1/(1-0.1666666667)

##plot 
plot.loh.tcga.bridges.genes.all.tumor.types.freq <- ggplot(stats.loh.tcga.all.tumor.types.AF,aes(x=bridges_genes,y=round(freq_AF*100,digits = 1), fill=Type))+ 
  geom_bar(stat="identity")+scale_fill_hue()+
  #geom_point()+scale_fill_hue()+ 
  theme(axis.text.x = element_text(size=8,angle = 90, vjust = 1,hjust=1))+
  scale_y_continuous(labels = function(x) paste0(x, "%"),breaks=c(0,20,40,60,80,100), limits=c(0,100))+
  facet_wrap(~Type, scale="free")
ggsave(filename = paste0("figures/plot.loh.tcga.bridges.genes.all.tumor.types.freq.pdf"), plot = plot.loh.tcga.bridges.genes.all.tumor.types.freq, width = 48, height = 28)

####plot EI
plot.EI.loh.tcga.bridges.genes.all.tumor.types <- ggplot(stats.loh.tcga.all.tumor.types.AF,aes(x=bridges_genes,y=EI, fill=Type))+ 
  #geom_bar(stat="identity")+scale_fill_hue()+ 
  geom_point(size = 4)+scale_fill_hue(l = 40, c = 30)+ 
  theme(axis.text.x = element_text(size=12,angle = 90, vjust = 1,hjust=1),axis.text.y = element_text(size=12,angle = 90, vjust = 1,hjust=1))+
  scale_y_continuous(labels = function(x) paste0(x),limits=c(1,2)) +
facet_wrap(~Type, scale="free") #scale_y_continuous(labels = function(x) paste0(x, "%"),breaks=c(0,20,40,60,80,100), limits=c(0,100))

ggsave(filename = paste0("plot.EI.loh.tcga.bridges.genes.all.tumor.types_updated.pdf"), plot = plot.EI.loh.tcga.bridges.genes.all.tumor.types, width = 48, height = 28)

plot.EI.loh.tcga.bridges.genes.BRCA <- ggplot(stats.loh.tcga.all.tumor.types.AF %>% filter (Type=="BRCA"),aes(x=bridges_genes,y=EI, fill=Type))+ 
  geom_point(size = 8)+scale_fill_hue()+ theme(axis.text.x = element_text(size=22,angle = 90, vjust = 1,hjust=1),
                                               axis.text.y = element_text(size=22,angle = 90, vjust = 1,hjust=1),
                                               axis.title=element_text(size=22)) +scale_y_continuous(labels = function(x) paste0(x),limits=c(1,2))
ggsave(filename = paste0("plot.EI.loh.tcga.bridges.genes.BRCA_updated_y-axis.pdf"), plot = plot.EI.loh.tcga.bridges.genes.BRCA, width = 48, height = 28)


