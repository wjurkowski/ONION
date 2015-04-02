
## @knitr include=FALSE
library(ReactomePA)
library(clusterProfiler)
library(biomaRt)
library(knitr)
library(org.Dr.eg.db)
opts_chunk$set(tidy=TRUE,tidy.opts=list(keep.blank.line=FALSE, width.cutoff=50),out.truncate=80,out.lines=6,cache=TRUE,dev='pdf',include=TRUE,fig.width=6,fig.height=6,resolution=150)
## @knitr options,results='hide',echo=FALSE
options(digits=3, width=80, prompt=" ", continue=" ")
args <- commandArgs(trailingOnly = TRUE)

## @knitr load sample data
#geneList <- as.data.frame(args[1])
gall <- read.table("~/Documents/Praca/Projects//Done/Simulife/nutrimouse/results/PLS-all.txt", header=FALSE)
g1 <- read.table("~/Documents/Praca/Projects//Done/Simulife/nutrimouse/results/PLS-1.txt", header=FALSE)
g2 <- read.table("~/Documents/Praca/Projects//Done/Simulife/nutrimouse/results/PLS-2.txt", header=FALSE)
g3 <- read.table("~/Documents/Praca/Projects//Done/Simulife/nutrimouse/results/PLS-3.txt", header=FALSE)
g4 <- read.table("~/Documents/Praca/Projects//Done/Simulife/nutrimouse/results/PLS-4.txt", header=FALSE)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
getEntrez <-function(ge){
  res <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), filters = "hgnc_symbol", values = ge, mart = mart)
  chr<-as.character(res$entrezgene)
  return(chr)
}
eall<-getEntrez(g4)
path="~/Documents/Praca/Projects/Done/Simulife/nutrimouse/results/"
write.table(eall,paste(path,"list4.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

require(ReactomePA)
xall <- enrichKEGG(gene=eall,pvalueCutoff=0.5, readable=T)
head(summary(xall))
# define biomart object
mart <- useMart(biomart = "REACTOME", dataset="reaction")
# query biomart
results <- getBM(attributes = c("stableidentifier_identifier","referencedatabase_chebi", "referencedatabase_uniprot", "referencedatabase_ensembl"), filters = c("referencemolecule_chebi_id_list","species_selection"), values = list(chebi,"Homo sapiens"), mart = mart)
av <- results[, "referencedatabase_ensembl", drop=FALSE]

## @knitr barplot, fig.cap="barplot of Reactome Pathway enrichment result.", fig.align="center", fig.height=5, fig.width=8, out.width="0.8\\textwidth", fig.pos="h"
barplot(x, showCategory=8)


## @knitr cnetplot, fig.cap="cnetplot of Reactome Pathway enrichment result.", fig.align="center", fig.height=16, fig.width=16, out.width="0.9\\textwidth", fig.pos="h"
cnetplot(x, categorySize="pvalue", foldChange=geneList)


## @knitr clusterProfiler, fig.cap="ReactomePA with clusterProfiler.", fig.align="center", fig.height=6, fig.width=11, out.width="0.9\\textwidth", fig.pos="h"
require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")
plot(res)


## @knitr GSEA analysis
y <- gseAnalyzer(geneList, nPerm=100, 
                  minGSSize=120, pvalueCutoff=0.05, 
                  pAdjustMethod="BH", verbose=FALSE)
res <- summary(y)
head(res)


## @knitr gseaplot, fig.cap="plotting gsea result", fig.align="center", fig.height=6, fig.width=6, out.width="0.6\\textwidth", fig.pos="h"
topID <- res[1,1]
topID
plot(y, geneSetID = topID)


## @knitr viewPathway, fig.cap="Reactome Pathway visualization.", fig.align="center", fig.height=16, fig.width=16, out.width="0.9\\textwidth", fig.pos="h"
viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)


## @knitr sessInfo, results='asis', echo=FALSE
toLatex(sessionInfo())


