exprSet <- read.table("genes.TMM.EXPR.matrix")
exprSet<- exprSet[rowSums(exprSet!= 0) > 0, ]
if(F){
#Convert ENSG to a gene symbol,Method 1
# Install and load the biomaRt package
library(biomaRt)
rownames(exprSet)
convert_ensg_to_gene_symbol <- function(ensg_list) {
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = ensg_list, mart = ensembl)
  gene_info <- gene_info[which(gene_info$external_gene_name != ""),]
  return(gene_info)
  ensg_list <-rownames(exprSet)
  gene_symbols <- convert_ensg_to_gene_symbol(ensg_list)
  print(gene_symbols)
}

#############################Method 2
library(tidyverse)
library(clusterProfiler)

count <-rownames(exprSet)
name <- bitr(count,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
}
exprSet=as.matrix(exprSet)
ensem2gene=read.csv(file="ensem2gene.csv")

length(unique(ensem2gene[ ,2]))
###Probes without gene annotation were removed
table(rownames(exprSet) %in% ensem2gene[ ,1])
#exprSet <- exprSet[ rownames(exprSet) %in% ensem2gene[ ,1], ]
#dim(exprSet)
ID2gene <- ensem2gene[ match(rownames(exprSet), ensem2gene[ ,1] ), ]
dim(ID2gene)
length(unique(ID2gene$gene_name))

rownames(exprSet) <- ID2gene[match(rownames(exprSet),ID2gene[,1]),2]
dim(exprSet)
length(unique(rownames(exprSet)))
########Check if a gene corresponds to more than one probe
tail( sort( table( ID2gene[ , 2 ] ) ), n = 12L )
##The expression data of the same gene is taken as the maximum value
MAX = by( exprSet, ID2gene[ , 2 ], function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX = as.character(MAX)
exprSet = exprSet[ rownames(exprSet) %in% MAX , ]
#########Removing duplicate genes
n_exprSet<-data.frame(exprSet,rownames(exprSet))
dim(n_exprSet)
if(F)################Determine whether there are annotated genes
{n_exprSet.na <- n_exprSet[is.na(n_exprSet$rownames.exprSet.)]
dim(n_exprSet.na)
exprset3=n_exprSet
failed <- exprset3[exprset3$Symbol=='',]
dim(failed)
}
microexprgenes.differ = n_exprSet
{
  microexprgenes.differ.dup <- microexprgenes.differ[duplicated(microexprgenes.differ$rownames.exprSet.),]
  dim(microexprgenes.differ.dup)
  microexprgenes.differ <- microexprgenes.differ[!duplicated(microexprgenes.differ$rownames.exprSet.),]
  dim(microexprgenes.differ)
}
exprSet1<-microexprgenes.differ
exprSet<-exprSet1[,-31]
dim(exprSet)

########The box plot detects the sample situation
{
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  boxplot(exprSet, col = cols,main="expression value",las=2)
}

###############################normalize
library(limma)
exprSet=normalizeBetweenArrays(exprSet)

ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0),))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
#ex[which(ex == 0)] <- 0.001
ex <- ex+ 1
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
{
  par(cex = 0.7)
  n.sample=ncol(exprSet1)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  boxplot(exprSet, col = cols,main="expression value",las=2)
}


write.csv(exprSet,file="exprSet.csv")
