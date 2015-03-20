options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)
library(biomaRt)
#library(data.table)

#read data
genes <- read.table(args[1], header=FALSE)
#w <- read.table(args[2], header=TRUE)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#iterate vector of ENSEMBL IDs
vec <- vector()
for(i in genes){
        # query biomart
        resu <- getBM(attributes = c("ensembl_peptide_id", "refseq_mrna","hgnc_symbol","entrezgene"), filters = "ensembl_peptide_id", values = i, mart = mart)
        # vec <- c(vec, results[,2])
}
#res <- unique(vec[vec != ""])
results <- as.data.frame(sapply(resu,gsub,pattern="^$",replacement="NA"))
if (args[2] == 'hgnc'){
  res<-unique(results$hgnc_symbol)
} else {
   res<-unique(results$refseq_mrna)
}

#d1 <- as.data.table(results)
#res <- d1[,list(refseq_mrna = paste(unique(refseq_mrna),collapse = ',')),by = 'ensembl_peptide_id']
#d2 <- as.data.table(w)
#setnames(res,"ensembl_peptide_id","node1_external_id")
#d3 <- merge(d2,res,by='node1_external_id', all.x=TRUE)
#setnames(res,"node1_external_id","node2_external_id")
#d4 <- merge(d3,res,by='node2_external_id', all.x=TRUE)

write.table(res,paste(args[1], args[2], sep='.'),row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
#write.table(d4,paste(args[2], '.mapped.txt', sep=''),row.names=FALSE,quote=FALSE,sep="\t")

