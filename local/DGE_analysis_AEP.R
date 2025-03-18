##### Hydra vulgaris re-analysis from dataset in Cazet et al. eLife, 2021 #####

#load annotations
Annotations <- read.csv("resources/dovetail_SP_annot.csv", stringsAsFactors = F, row.names = 1)

#load read counts (genome)

counts <- read.table("resources/AEP_mat.txt")

#rename sample names in columns and reorder
colnames(counts) <- c("F0R1", "F0R2", "F0R3",
                      "H0R1", "H0R2", "H0R3",
                      "F12R1","F12R2","F12R3",
                      "H12R1", "H12R2", "H12R3",
                      "F3R1", "F3R2", "F3R3",
                      "H3R1","H3R2","H3R3",
                      "F8R1","F8R2","F8R3",
                      "H8R1","H8R2","H8R3",
                      "iF0R1", "iF0R2", "iF0R3",
                      "iH0R1", "iH0R2", "iH0R3",
                      "iF12R1","iF12R2","iF12R3",
                      "iH12R1", "iH12R2", "iH12R3",
                      "iF8R1","iF8R2","iF8R3",
                      "iH8R1","iH8R2","iH8R3")

counts <- counts[,!grepl("i",colnames(counts))]

#add row with gene IDs
counts$ID <- rownames(counts)

#define treatment groups
TR <- factor(c(rep("F0", 3), rep("H0", 3),
               rep("F12", 3), rep("H12", 3),
               rep("F3", 3), rep("H3", 3),
               rep("F8", 3), rep("H8", 3)))

#initiate DGElist object
analysisList <- DGEList(counts = counts[,1:(ncol(counts)-1)], group = TR, genes = counts$ID)

#exclude lowly expressed genes from the analysis 
#(at least three samples have at least two counts per million)

keep <- rowSums(cpm(analysisList)>=2) >= 3
analysisList <- analysisList[keep, , keep.lib.sizes=F]

#calculate normalization factors for each sample
analysisList <- calcNormFactors(analysisList)

#define treatment groups
design <- model.matrix(~0+TR, data = analysisList$samples)

colnames(design) <- gsub("TR","",colnames(design))

#calculate differential gene expression
analysisList <- estimateGLMRobustDisp(analysisList,design)

fit <- glmFit(normdge, design, robust = T)


#define pairwise comparisons to be used
my.contrasts <- makeContrasts(
  
  H3vH0 = `H3` - `H0`,
  H8vH0 = `H8` - `H0`,
  H12vH0 = `H12` - `H0`,
  
  F3vF0 = `F3` - `F0`,
  F8vF0 = `F8` - `F0`,
  F12vF0 = `F12` - `F0`,
  
  F8vF3 = `F8` - `F3`,
  F12vF3 = `F12` - `F3`,
  
  H8vH3 = `H8` - `H3`,
  H12vH3 = `F8` - `F3`,
  
  H12vH8 = `H12` - `H8`,
  F12vF8 = `F12` - `F8`,
  
  HR3vFR3 = (`H3` - `H0`) - (`F3` - `F0`),
  HR8vFR8c0 = (`H8` - `H0`) - (`F8` - `F0`),
  HR8vFR8c3 = (`H8` - `H3`) - (`F8` - `F3`),
  
  HR8vFR12c0 = (`H8` - `H0`) - (`F12` - `F0`),
  HR8vFR12c3 = (`H8` - `H3`) - (`F12` - `F3`),
  
  HR12vFR12c0 = (`H12` - `H0`) - (`F12` - `F0`),
  HR12vFR12c3 = (`H12` - `H3`) - (`F12` - `F3`),
  
  levels=design)

#Define thresholds for glmTreat

fc = 0
fdr = 5e-2

####Perform Pairwise Comparisons####

genResTables <- function(x,fdr,fc) {
  
  y <- glmTreat(fit, contrast = my.contrasts[,x], lfc = log2(1.2), null = "interval")
  y <- topTags(y, n=Inf)
  #extract results table
  z <-  y$table
  z$ID <- rownames(z)
  z <- merge(z,Annotations, by= "ID", all.x = T)
  
  #significance filters
  z.DG <- z[z$FDR <= fdr,]
  z.DG.up <- z.DG[z.DG$logFC > fc,]
  z.DG.down <- z.DG[z.DG$logFC < -fc,]
  
  #return the four results objects
  return(list(z,z.DG,z.DG.up,z.DG.down))
  #return(list(z))
}

#generate all results objects
#systematize object naming based on contrast names

for (i in 1:length(colnames(my.contrasts))) {
  
  resList <- genResTables(colnames(my.contrasts)[i],fdr,fc)
  
  assign(colnames(my.contrasts)[i], resList[[1]])
  assign(paste0(colnames(my.contrasts)[i],".DG"), resList[[2]])
  assign(paste0(colnames(my.contrasts)[i],".DG.up"), resList[[3]])
  assign(paste0(colnames(my.contrasts)[i],".DG.down"), resList[[4]])
  
}

####Save the results####

#We also want to save the normalized counts (for plotting purposes)
normCountsAEP <- as.data.frame(cpm(analysisList, normalized.lib.sizes=T, log = T))
normCountsAEP$ID <- rownames(normCountsAEP)
normCountsAEP <- merge(normCountsAEP,Annotations, by="ID", all.x = T)
write.csv(normCountsAEP,file = "DGE_tables_AEP/normCounts.AEP.csv")

#save results dataframes, the DGE object, Annotations, and the normalized counts

#save results from all contrasts
saveThese <- ls(pattern = paste(colnames(my.contrasts),collapse = "|"))

# save gene lists in a csv (more git friendly format)
lapply(saveThese, function(x) write.csv(eval(parse(text=x)), file = paste0("DGE_tables_AEP/",x,".csv")))

#save contrasts table to help make sense of file names
write.csv(my.contrasts, file = "DGE_tables_AEP/contrasts.csv")

saveThese <- c(saveThese,"analysisList","normCountsAEP","Annotations")

save(file = "DGE_tables_AEP/RNA_DGE.RData", list = saveThese)
