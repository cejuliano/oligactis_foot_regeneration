### Install matrixStats package
install.packages("matrixStats")
library(matrixStats)
install.packages("Hmisc")
library("Hmisc")	

flattenCorrMatrix <- function(cormat, pmat) {
	ut <- upper.tri(cormat)
	data.frame(row=rownames(cormat)[row(cormat)[ut]],
	column=rownames(cormat)[col(cormat)[ut]], cor=(cormat)[ut], p=pmat[ut])
}

### For loop to calculate Co-expression Network for each gene expression profile.


# change this according to the input expression matrix.
# assuming the input files is in the "results" directory.
Input.oli = "resources/fpkm.OLI.csv"

Input.vul = "resources/fpkm.AEP.csv"


### Read input file 
GE.oli=read.table(Input.oli, sep = ",", header = T,stringsAsFactors= F, row.names=1)
keep <- rowSums(GE.oli>=2) >= 3
GE.oli <- GE.oli[keep, ]
GE.vul=read.table(Input.vul, sep = ",", header = T,stringsAsFactors= F, row.names=1)
keep <- rowSums(GE.vul>=2) >= 3
GE.vul <- GE.vul[keep, ]
GE.tmp.oli=GE.oli
GE.tmp.vul=GE.vul

### The input data can be filtered based on their expression patterns. 
### Remove genes with zero standard deviation.
### SD != 0
library(matrixStats)
GE.tmp.oli=as.matrix(GE.oli)
GE.tmp.vul=as.matrix(GE.vul)
GE.rowSds.oli = transform(GE.tmp.oli, SD=rowSds(GE.tmp.oli, na.rm=TRUE))
GE.rowSds.vul = transform(GE.tmp.vul, SD=rowSds(GE.tmp.vul, na.rm=TRUE))
GE.rowSds.oli = GE.rowSds.oli[which(GE.rowSds.oli [,40]!=0),]
GE.rowSds.vul = GE.rowSds.vul[which(GE.rowSds.vul [,25]!=0),]
GE.rowSds.oli = GE.rowSds.oli[order(GE.rowSds.oli [,40]),]
GE.rowSds.vul = GE.rowSds.vul[order(GE.rowSds.vul [,25]),]
GE.tmp.oli = GE.rowSds.oli[, (1:39)]
GE.tmp.vul = GE.rowSds.vul[, (1:24)]
		

### Keep only genes with variance >1
GE.tmp.oli=as.matrix(GE.tmp.oli)
GE.tmp.vul=as.matrix(GE.tmp.vul)
GE.rowVars.oli = transform(GE.tmp.oli, Vars=rowVars(GE.tmp.oli, na.rm=TRUE))
GE.rowVars.vul = transform(GE.tmp.vul, Vars=rowVars(GE.tmp.vul, na.rm = TRUE))
GE.rowVars.5.oli = GE.rowVars.oli[which(GE.rowVars.oli[,40] > 1),]
GE.rowVars.5.vul = GE.rowVars.vul[which(GE.rowVars.vul[,25] > 1),]
GE.rowVars.5.oli = GE.rowVars.5.oli[order(row.names(GE.rowVars.5.oli)), ]
GE.rowVars.5.vul = GE.rowVars.5.vul[order(row.names(GE.rowVars.5.vul)), ]
GE.tmp.oli=GE.rowVars.5.oli[, -40]
GE.tmp.vul=GE.rowVars.5.vul[, -25] 
	
### Calculating correlation between genes with corr & P-val
GE_gene_rcor.oli <- rcorr(as.matrix(t(GE.tmp.oli)))
GE_gene_rcor.vul <- rcorr(as.matrix(t(GE.tmp.vul)))

### Converting cor matrix to cor table by flattenCorrMatrix function
GE_gene_rcor_list.oli =flattenCorrMatrix(GE_gene_rcor.oli$r, GE_gene_rcor.oli$P)
GE_gene_rcor_list.vul =flattenCorrMatrix(GE_gene_rcor.vul$r, GE_gene_rcor.vul$P)


### Filtering rcor result
pvalCut=0.001
corCut=0.975

### p-val cutoff
GE_gene_rcor_list_tmp.oli = GE_gene_rcor_list.oli
GE_gene_rcor_list_tmp.vul = GE_gene_rcor_list.vul
GE_gene_rcor_list_filterP.oli= GE_gene_rcor_list_tmp.oli[which(GE_gene_rcor_list_tmp.oli$p<=pvalCut),]
GE_gene_rcor_list_filterP.vul= GE_gene_rcor_list_tmp.vul[which(GE_gene_rcor_list_tmp.vul$p<=pvalCut),]
GE_gene_rcor_list_tmp.oli = GE_gene_rcor_list_filterP.oli
GE_gene_rcor_list_tmp.vul = GE_gene_rcor_list_filterP.vul

### cor cutoff
GE_gene_rcor_list_filterCor.oli = GE_gene_rcor_list_tmp.oli[which(GE_gene_rcor_list_tmp.oli$cor>=corCut),]
GE_gene_rcor_list_filterCor.vul = GE_gene_rcor_list_tmp.vul[which(GE_gene_rcor_list_tmp.vul$cor>=corCut),]
GE_gene_rcor_list_tmp.oli= GE_gene_rcor_list_filterCor.oli
GE_gene_rcor_list_tmp.vul= GE_gene_rcor_list_filterCor.vul
GE_gene_rcor_list_Out.oli = cbind(GE_gene_rcor_list_tmp.oli[,1:2], cor=round(GE_gene_rcor_list_tmp.oli[,3],3), pval=signif(GE_gene_rcor_list_tmp.oli[,4],3))
GE_gene_rcor_list_Out.vul = cbind(GE_gene_rcor_list_tmp.vul[,1:2], cor=round(GE_gene_rcor_list_tmp.vul[,3],3), pval=signif(GE_gene_rcor_list_tmp.vul[,4],3))	
	
### Writing filtered rcor result into csv file
Output.oli=paste0("CoexpNetwork_", Input.oli)
Output.vul=paste0("CoexpNetwork_", Input.vul)
write.csv(GE_gene_rcor_list_Out.oli, file=Output.oli)
write.csv(GE_gene_rcor_list_Out.vul, file=Output.vul)

### Reading table with reciprocal blast hits
RBH = "HOLIAEP.RBH.txt"
rbh.table <- read.table(RBH, sep = ",", header = F, stringsAsFactors = F)
rbh.ortho <- rbh.table[,1:2]
### Filter rbh table by AEP genes for which we have expression data only
tmp <- data.frame(rbh.ortho)
colnames(tmp) <- c("AEP","OLI")
tmp2 <- data.frame(rownames(GE.vul))
colnames(tmp2) <- "AEP"  
tmp3 <- merge(tmp,tmp2, by="AEP")
rbh.ortho <- tmp3
  
### Loading network data
OLI_GE_gene_rcor_list = read.table("CoexpNetwork_fpkm.OLI.csv", sep= ",", header=TRUE,row.names=1,as.is=T)
AEP_GE_gene_rcor_list = read.table("CoexpNetwork_fpkm.AEP.csv", sep=",", header=TRUE,row.names=1,as.is=T)

### Co-expression input files for OrthoClust R
OLI_edgelist = as.matrix(OLI_GE_gene_rcor_list[,1:2])
#colnames(OLI_edgelist) = C()
AEP_edgelist = as.matrix(AEP_GE_gene_rcor_list[,1:2])
  
KAPPA=2
ortho.results=OrthoClust2(Eg1 = AEP_edgelist, Eg2 = OLI_edgelist, list_orthologs = rbh.ortho, kappa = KAPPA)

### OrthoClust result into csv file
AEP_modules=as.matrix(ortho.results[[1]])
OLI_modules=as.matrix(ortho.results[[2]])
AEP_modules=as.data.frame(cbind(AEP_modules, "AEP", rownames(AEP_modules)),stringsAsFactors=FALSE)
OLI_modules=as.data.frame(cbind(OLI_modules, "OLI", rownames(OLI_modules)),stringsAsFactors=FALSE)
Merged_modules=rbind(AEP_modules, OLI_modules)
colnames(Merged_modules) = c("ModuleID", "Species", "GeneID")
Merged_modules[, 1] <- as.numeric(Merged_modules[, 1])
Merged_modules = Merged_modules[order(Merged_modules$ModuleID),]
rownames(Merged_modules) = NULL


### Organizing module info into a module number table 
NumModule =as.numeric(rownames(as.matrix(tail(table(ortho.results$species1_clust), n=1))));
ortho_module_list = matrix(0, nrow = NumModule, ncol = 4)
ortho_module_list = cbind(as.matrix(1:NumModule), ortho_module_list) 
ortho_tb1=table(ortho.results$species1_clust)
ortho_tb2=table(ortho.results$species2_clust)
ortho_tb1_list = cbind(Module=as.numeric(rownames(as.matrix(ortho_tb1))), NoGenes=as.numeric(as.matrix(ortho_tb1)))
ortho_tb2_list = cbind(Module=as.numeric(rownames(as.matrix(ortho_tb2))), NoGenes=as.numeric(as.matrix(ortho_tb2)))

### Filling in the module number table 
for (i in 1:NumModule)
{
	idmx1 = match(ortho_module_list[i,1], ortho_tb1_list[,1])
	modNum1= ortho_tb1_list[idmx1,1]
	numGen1= ortho_tb1_list[idmx1,2]
	ortho_module_list[i,2]=modNum1
	ortho_module_list[i,3]=numGen1
	
	idmx2 = match(ortho_module_list[i,1], ortho_tb2_list[,1])
	modNum2= ortho_tb2_list[idmx2,1]
	numGen2= ortho_tb2_list[idmx2,2]
	ortho_module_list[i,4]=modNum2
	ortho_module_list[i,5]=numGen2
}
ortho_module_list=ortho_module_list[, c(1,3,5)]
ortho_module_list[is.na(ortho_module_list)] <- 0

### Calculating percentages for the module number table
#install.packages("matrixStats")
library(matrixStats)
tmp = ortho_module_list
ortho_module_list=cbind(ortho_module_list,sum=ortho_module_list[,2]+ortho_module_list[,3])
ortho_module_list=cbind(ortho_module_list,ratio1=round(ortho_module_list[,2]/ortho_module_list[,4],3), ratio2=round(ortho_module_list[,3]/ortho_module_list[,4],3))
colnames(ortho_module_list)=c("ModuleID","NumOfGenesFromSpe1","NumOfGenesFromSpe2","TotalNumOfGenesFromModule","RatioOfGenesFromSpe1","RatioOfGenesFromSpe2")
ortho_module_list_sumOrdered <- ortho_module_list[order(ortho_module_list[,4],decreasing = TRUE),] 


### Saving OrthoClust result into csv
OrthoOutFile=paste0("Orthoclust_Results.csv")
write.csv(Merged_modules , OrthoOutFile)

### Saving a table of OrthoClust result into csv file
TableFile=gsub(".csv", "_Summary.csv", OrthoOutFile)
write.csv(ortho_module_list_sumOrdered , TableFile)

### Saving OrthoClust results into RData
OrthoRDataFile=gsub(".csv", ".RData", OrthoOutFile)
save(ortho_module_list_sumOrdered, ortho.results, AEP_edgelist, OLI_edgelist, rbh.ortho, KAPPA, file = OrthoRDataFile)
	
### Saving co-expression network into RData
RDataFile=paste0("CoexpNetwork_", gsub(".csv", "", Input.oli),".RData")
print("RDataFile name:"); print(RDataFile)
save(GE.oli, GE.tmp.oli, GE_gene_rcor_list_Out.oli, file = RDataFile)
save(GE.vul, GE.tmp.vul, GE_gene_rcor_list_Out.vul, file = RDataFile)