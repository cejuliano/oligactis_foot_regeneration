if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("maSigPro")

library(maSigPro)

Input <- "mat.ALP.OLI.csv" #Matrix containing ID and cpm per timepoint 
data.mat <- read.table(Input,header = TRUE, sep = ",",row.names = 1)
data.mat <- 2^data.mat #transform logarithmic counts to linear counts per million
Input <- "design.ALP.OLI.csv" #File with treatment and replicates information
design.mat <- read.table(Input, header = TRUE, sep = ",", row.names = 1)

design <- make.design.matrix(design.mat, degree = 4) #change degrees according to number of time points

#FILTER STEP: filter out influential genes from matrix SKIP in the first round
data.mat$ID <- rownames(data.mat)
data.mat <- data.mat[!(data.mat$ID %in% search.vector),]
data.mat$ID <- NULL

#Find genes that change their expression level significantly across the time course
fit <- p.vector(data.mat, design, counts = TRUE, theta = 10, Q = 0.05, MT.adjust = "BH")#Q=0.001
fit$i

#Find significant variables for each gene
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05) #0.01

search.vector <- colnames(tstep$influ.info) #Go back to FILTER STEP

#Obtain list of significant genes
sigs <- get.siggenes(tstep, rsq = 0.80, vars = "groups") #R-squared may be adjusted to modulate pattern consistency vs number of signifcant genes
sigs$sig.genes$ALPvsControl$g

#
p <- see.genes(sigs$sig.genes$ALPvsControl, show.fit = F, dis =design$dis,
               cluster.method="hclust" ,cluster.data = 1, k.mclust = TRUE)

clusters <- as.data.frame(names(p$cut))
clusters$group <- as.numeric(p$cut)
colnames(clusters) <- c("ID","cluster.mSP")

annotations <- read.csv("annotations.csv", header = TRUE, sep = ",")

cluster.list <- merge(clusters, annotations, by = "ID")
write.csv(cluster.list, file = "cluster.list.csv")