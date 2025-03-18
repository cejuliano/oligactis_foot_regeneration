###Install Packages required#####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

install.packages("rstudioapi")

install.packages("ggfortify")

install.packages('gsubfn')

library(sva)
library(edgeR)
library(limma)
library(rstudioapi)
library(ggfortify)
library(gsubfn)

#Open raw counts matrix

Input = "resources/OLI_mat.txt"
counts <- read.table(Input, header=TRUE, row.names = 1, sep="\t", comment.char="")
rowid <- row.names(counts)
counts.names <- colnames(counts)

#Create a vector to indicate which samples belong to which sequencing batch

batch <- c(1,1,1,1,2,2,2,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,
           2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,2,2)

#Create a vector to indicate the experimental group (time point) to which each sample belongs

group <- c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,
           11,11,11,12,12,12,13,13,13,14,14,14,14,15,15,15,15,15,16,16,16,16,16)

#Reformating of the data frame to work with ComBat_seq

counts<- data.matrix(counts, rownames.force = NA)

#Batch normalization 

counts <- ComBat_seq(counts, batch = batch, group = group)

#Restructure counts matrix to match sample (columns) and gene (rows) number 

counts <- matrix(data = counts, nrow =23396, ncol = 61)

colnames(counts) = counts.names

rownames(counts) <- rowid

counts <- as.data.frame(counts)

#Add IDs as column

counts = counts[,!grepl("i",colnames(counts))]

counts$ID <- rownames(counts)

#define treatment groups

factor1 = substr(colnames(counts[1:length(counts)-1]), 1, 1)
factor2 = substr(colnames(counts[1:length(counts)-1]), 4, 4)

grp = interaction(factor1,factor2)

#Extract group names leaving the replicate number out

grp = sub("..$", "", colnames(counts[1:length(counts)-1]))

#Create the Differential Gene Expression object

dge <- DGEList(counts = counts[,1:(ncol(counts)-1)], group = grp, genes = counts$ID)

#exclude lowly expressed genes from the analysis 
#(at least six samples have at least two counts per million)

keep <- rowSums(cpm(dge)>=2) >= 3

dge <- dge[keep, , keep.lib.sizes=F]

#calculate normalization factors for each sample

normdge <- calcNormFactors(dge)
normdge = estimateCommonDisp(normdge)
normdge$common.dispersion

#Define groups for DGE analysis

design <- model.matrix(~0+grp, data = normdge$samples) 

colnames(design) <- gsub("grp","",colnames(design))

#calculate differential gene expression

normdge <- estimateGLMRobustDisp(normdge,design)

fit <- glmFit(normdge, design, robust = T)

#define pairwise comparisons to be used
my.contrasts <- makeContrasts(
  
  F3vF0 = `FR_3H` - `FR_0H`,
  F12vF0 = `FR_12H` - `FR_0H`,
  F24vF0 = `FR_24H` - `FR_0H`,
  F48vF0 = `FR_48H` - `FR_0H`,
  F48Vfoot= `FR_48H` - `FOOT_O`,
  
  H3vH0 = `HR_3H` - `HR_0H`,
  H12vH0 = `HR_12H` - `HR_0H`,
  H24vH0 = `HR_24H` - `HR_0H`,
  H48vH0 = `HR_48H` - `HR_0H`,
  H48vhead = `HR_48H` - `HEAD_O`,
  
  ALP3vALP0 = `AL_3H` - `FR_0H`,
  ALP12vALP0 = `AL_12H` - `FR_0H`,
  ALP24vALP0 = `AL_24H` - `FR_0H`,
  ALP48vALP0 = `AL_48H` - `FR_0H`,
  ALP48vfoot = `AL_48H` - `FOOT_O`, 
  
  HR3vFR3 =  (`HR_3H` - `HR_0H`) - (`FR_3H` - `FR_0H`),
  HR12vFR12c0 = (`HR_12H` - `HR_0H`) - (`FR_12H` - `FR_0H`),
  HR24vFR24c0 = (`HR_24H` - `HR_0H`) - (`FR_24H` - `FR_0H`),
  HR48vFR48c0 = (`HR_48H` - `HR_0H`) - (`FR_48H` - `FR_0H`),
  
  HR3vAR3 =  (`HR_3H` - `HR_0H`) - (`AL_3H` - `FR_0H`),
  HR12vAR12c0 = (`HR_12H` - `HR_0H`) - (`AL_12H` - `FR_0H`),
  HR24vAR24c0 = (`HR_24H` - `HR_0H`) - (`AL_24H` - `FR_0H`),
  HR48vAR48c0 = (`HR_48H` - `HR_0H`) - (`AL_48H` - `FR_0H`),
  
  FR3vAR3 =  (`FR_3H` - `FR_0H`) - (`AL_3H` - `FR_0H`),
  FR12vAR12c0 = (`FR_12H` - `FR_0H`) - (`AL_12H` - `FR_0H`),
  FR24vAR24c0 = (`FR_24H` - `FR_0H`) - (`AL_24H` - `FR_0H`),
  FR48vAR48c0 = (`FR_48H` - `FR_0H`) - (`AL_48H` - `FR_0H`),
  
  FootvHead = `FOOT_O` - `HEAD_O`,
  FootvF0 = `FOOT_O` - `FR_0H`,
  FootvH0 = `FOOT_O` - `HR_0H`,
  HeadvF0 = `HEAD_O` - `FR_0H`,
  HeadvH0 = `HEAD_O` - `HR_0H`,
  
  levels=design)

# Define thresholds for glmTreat function

fc = 0
fdr = 1e-2

genResTables <- function(x,fdr,fc) {
  #pull full results
  y <- glmTreat(fit, contrast = my.contrasts[,x], lfc = log2(1.2), null = "interval")
  y <- topTags(y, n=Inf)
  #extract results table
  z <-  y$table
  z$ID <- rownames(z)
  z <- merge(z,annotations, by= "ID", all.x = T)
  #significance filters
  z.DG <- z[z$FDR <= fdr,]
  z.DG.up <- z.DG[z.DG$logFC > fc,]
  z.DG.down <- z.DG[z.DG$logFC < -fc,]
  
  #return the four results objects
  
  return(list(z,z.DG,z.DG.up,z.DG.down))
}

Input = "resources/OLI.AEP.annotation.csv" #Read annotations to add to results tables
annotations <- read.csv(Input, header = TRUE, sep = ",", stringsAsFactors = F, row.names = NULL)

#generate all results objects
#systematize object naming based on contrast names

for (i in 1:length(colnames(my.contrasts))) {
  
  resList <- genResTables(colnames(my.contrasts)[i],fdr,fc)
  
  assign(colnames(my.contrasts)[i], resList[[1]])
  assign(paste0(colnames(my.contrasts)[i],".DG"), resList[[2]])
  assign(paste0(colnames(my.contrasts)[i],".DG.up"), resList[[3]])
  assign(paste0(colnames(my.contrasts)[i],".DG.down"), resList[[4]])
}

#Save the normalized counts (for plotting purposes) and add annotations

normalizedCounts <- as.data.frame(cpm(normdge, normalized.lib.sizes=T, log = T))
normalizedCounts$ID <- rownames(normalizedCounts)
normCounts <- merge(normalizedCounts, annotations, by="ID", all.x = T)

write.csv(normCounts, file = "DGE_tables_OLI/normCounts.OLI.csv")

#save results from all contrasts

saveThese <- ls(pattern = paste(colnames(my.contrasts),collapse = "|"))

# save gene lists in a csv (more git friendly format)

lapply(saveThese, function(x) write.csv(eval(parse(text=x)), file = paste0("DGE_tables_OLI/","",x,".csv")))

#save contrasts table to help make sense of file names

write.csv(my.contrasts, file = "DGE_tables_OLI/contrasts.csv")

saveThese <- c(saveThese,"normdge","normalizedCounts")

save(file = "DGE_tables_OLI/RNA_DGE.RData", list = saveThese)

