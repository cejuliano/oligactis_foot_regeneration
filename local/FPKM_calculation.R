install.packages("DGEobj.utils")
library(DGEobj.utils)

Input = "resources/dovetail.length.csv"
lengths.vul=read.table(Input, sep = ",", header = T,stringsAsFactors= F, row.names=NULL)

Input = "resources/normCounts.AEP.csv"
cpm.vul=read.table(Input, sep = ",", header = T,stringsAsFactors= F, row.names=NULL)

lengths.vul = merge(lengths.vul, cpm.vul, by = "ID", all.x = FALSE)

lengths.vul = as.data.frame(cbind(lengths.vul$ID,lengths.vul$length))
colnames(lengths.vul) <- c("ID","length")

#Transform counts back from log scale
tmp <- cpm.vul[,3:26]
tmp <- apply(tmp, 2, function(x) {
  2^x
})
#Obtain gene.Length vector
tmp2 <- sapply(lengths.vul$length,as.numeric)

fpkm.vul = convertCounts(tmp, unit = "fpkm",geneLength = tmp2)
row.names(fpkm.vul) = cpm.vul$ID
write.csv(fpkm.vul, file = "resources/fpkm.AEP.csv")
