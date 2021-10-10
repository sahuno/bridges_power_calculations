### compute loh in TCGA breast cacner copy number 
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


#get breast cancer sample names, from aneuploidy data
###subset breast cancer from tcga
df.tcga.aneuploidy <- fread("/Users/samuelahuno/Downloads/PANCAN_ArmCallsAndAneuploidyScore_092817.txt", sep = "\t", header= TRUE)
df.tcga.aneuploidy.brca <- df.tcga.aneuploidy %>% filter(Type=="BRCA")
head(df.tcga.aneuploidy.brca)
dim(df.tcga.aneuploidy.brca)

##is there NAs
is.na(df.tcga.aneuploidy.brca$Sample)
anyNA(df.tcga.aneuploidy.brca$Sample)

##check distinct samples
unique.brca <- distinct(df.tcga.aneuploidy.brca,Sample)
df.tcga.aneuploidy %>% filter(Sample=="TCGA-OR-A5J1-01")

merged.tcga.mastercalls.abs.BRCA.ID <- merge(df.tcga.mastercalls.abs, unique.brca,by.x="Sample",by.y="Sample",no.dups=TRUE)
dim(merged.tcga.mastercalls.abs.BRCA.ID)

df.tcga.mastercalls.abs[match(df.tcga.mastercalls.abs$Sample,unique.brca$Sample),]

merged.tcga.mastercalls.abs.BRCA.ID %>% distinct_at(vars(Sample),.keep_all=T)  %>% dplyr::summarise(n = n())
##subset only loh
lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- merged.tcga.mastercalls.abs.BRCA.ID %>% filter(LOH == 1)
distinct(lohmerged.df.tcga.mastercalls.abs.BRCA.ID,Sample)


#sanity check, see if there are samples, i missed
x=c("a","b","c")
y=c("b","a")
setdiff(x,y)
setdiff(df.tcga.aneuploidy.brca$Sample,merged.tcga.mastercalls.abs.BRCA.ID$Sample)
unique(df.tcga.aneuploidy.brca$Sample[! df.tcga.aneuploidy.brca$Sample %in% merged.tcga.mastercalls.abs.BRCA.ID$Sample])

### read breast cacner genes from bridges study
source("/Users/samuelahuno/Polak_lab10082019/bridges_breast/scripts/getCytobandsHg19_bridges_genes.R")
cordinates.bridges.genes.hg19

##find overlaps
gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- makeGRangesFromDataFrame(lohmerged.df.tcga.mastercalls.abs.BRCA.ID,start.field=  "Start", end.field= "End", seqnames.field="Chromosome", keep.extra.columns = T)
ov.gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- findOverlaps(query=cordinates.bridges.genes.hg19, subject=gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID, type = "within") ##complete overlaps

####append, gene names to columns
gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID$bridges_genes <- NA
mcols(gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID)$bridges_genes[subjectHits(ov.gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID)] <- mcols(cordinates.bridges.genes.hg19)$hgnc_symbol[queryHits(ov.gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID)]


##convert to data frame
df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- as.data.frame(gr.lohmerged.df.tcga.mastercalls.abs.BRCA.ID)
df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% drop_na(bridges_genes)
dim(df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID)


###might not need this line
df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID <- df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% distinct_at(vars(Sample),.keep_all=T) %>% group_by(bridges_genes) 
df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% filter()

##grab from pp scipts####
stats.loh.tcga.AF <- df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% group_by(bridges_genes) %>% dplyr::summarise(n = n()) %>% mutate(freq = n /1048,freq_AF = (n /1048)*0.5)
#view(stats.loh.tcga.AF)

bridges.genes.freq <- read.table(file="/Users/samuelahuno/Polak_lab10082019/bridges_breast/data/bridges.genes.freq.tsv",header = T,sep = '\t')
stats.loh.tcga.AF.MSig <- merge(stats.loh.tcga.AF,bridges.genes.freq[,c(1,4,5,6)],by.x="bridges_genes",by.y="gene")

stats.loh.tcga.AF.MSig$freq_AF <- round(stats.loh.tcga.AF.MSig$freq_AF*100,digits = 1)

stats.loh.tcga.AF.MSig <- stats.loh.tcga.AF.MSig[,-c(2,3,6)]
names(stats.loh.tcga.AF.MSig) <- c("Genes","%_biallellic_inactivation","Mutational_signatures","Proposed_therapy_choice")


write.table(stats.loh.tcga.AF.MSig, "data/stats.loh.tcga.AF.MSig.txt",quote=FALSE, sep='\t',na = "",row.names = FALSE,col.names = TRUE)

#HRD <- c("BRCA1","BRCA2","BARD1","RAD51C","RAD51D","PALB2") 
#HRD <- c("MLH1","PMS2","MSH6","MS2") 
#sig
#MLH1, PMS2 MSH6 and MS2 are for HRD signauttre
#MUTYH - Signature 18
##scratch
bridges.genes.freq <- read.table(file="/Users/samuelahuno/Polak_lab10082019/Administravia/Nejm.genes.txt",header = F)
write.table(bridges.genes.freq, "data/bridges.genes.freq.tsv",quote=FALSE, sep='\t',na = "",row.names = FALSE,col.names = FALSE)
##############



###stats for plot
stats.loh.tcga <- df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% group_by(bridges_genes) %>% dplyr::summarise(n = n()) %>% mutate(freq = n / sum(n),sumz=sum(n))
df.lohmerged.df.tcga.mastercalls.abs.BRCA.ID %>% group_by(bridges_genes) %>% dplyr::summarise(n = n())


view(stats.loh.tcga)
sum(stats.loh.tcga$n)

plot.loh.tcga.bridges.genes <- ggplot(stats.loh.tcga,aes(x=bridges_genes,y=freq))+geom_bar(stat="identity")+scale_fill_hue()
plot.loh.tcga.bridges.genes.numbers <- ggplot(stats.loh.tcga,aes(x=bridges_genes,y=n))+geom_bar(stat="identity")+scale_fill_hue()

ggsave(filename = paste0("figures/bar_plot_freq_loh_tcga_bridges_genes_v2.pdf"), plot = plot.loh.tcga.bridges.genes, width = 18, height = 9)
ggsave(filename = paste0("figures/bar_plot_numbers_loh_tcga_bridges_genes_v2.pdf"), plot = plot.loh.tcga.bridges.genes.numbers, width = 18, height = 9)


