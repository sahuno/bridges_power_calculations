## function to get gene cordinates, biomart
###date; Dec 09 2020 
library(biomaRt)
library(GenomeInfoDb)
library(GenomicFeatures)

####extracting cytopands and gene cordinates for comparison between tcga hg19 ensemsbl  ###############################
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# Only use standard human chromosomes
target.chroms <- c(1:22, "X", "Y", "M")

bridges_genes <- c("PIK3CA","MEN1","AKT1", "MLH1", "BABAM2","GEN1","RINT1","EPCAM","RECQL",
                   "CDH1","MRE11","NBN", "XRCC2","ABRAXAS1", "MUTYH", "MSH2", "FANCM", "RAD50", 
                   "BRIP1", "PMS2","FANCC", "STK11", "NF1", "RAD51D", "RAD51C", "MSH6", "BARD1", 
                   "ATM","PTEN", "CHEK2", "TP53", "PALB2", "BRCA2", "BRCA1")
length(bridges_genes)

##use ensembl ids  for genes
ensembl_ids <- c("ENSG00000158019","ENSG00000020922")

cordinates.bridges.genes.hg19 <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band","gene_biotype"),
                         filters = c("hgnc_symbol", "chromosome_name"),
                         values = list(hgnc_symbol=bridges_genes, chromosome_name=target.chroms),
                         mart = grch37)

cordinates.bridges.genes.hg19 <- makeGRangesFromDataFrame(cordinates.bridges.genes.hg19,start.field=  "start_position", end.field= "end_position", seqnames.field="chromosome_name", keep.extra.columns = T)


### what is th ewidth of genes
mcols(cordinates.bridges.genes.hg19)$width_genes <- width(cordinates.bridges.genes.hg19)

##sanity checks, di
length(mcols(cordinates.bridges.genes.hg19)$hgnc_symbol)
table(mcols(cordinates.bridges.genes.hg19)$hgnc_symbol)

found_genes <- mcols(cordinates.bridges.genes.hg19)$hgnc_symbol
setdiff(found_genes,bridges_genes)
unique(bridges_genes[! bridges_genes %in% found_genes])
#[1] "BABAM2"   "MRE11"    "ABRAXAS1" these genes not found, find other gene alias

#use hgnc numbers
# "BABAM2" = 
