options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)
library(biomaRt)

#read data
#chebi <- read.table(args[1], header=FALSE)
chebi <-as.data.frame(args[1])
chebi
# define biomart object
mart <- useMart(biomart = "REACTOME", dataset="reaction")
# query biomart
results <- getBM(attributes = c("stableidentifier_identifier","referencedatabase_chebi", "referencedatabase_uniprot", "referencedatabase_ensembl"), filters = c("referencemolecule_chebi_id_list","species_selection"), values = list(chebi,"Homo sapiens"), mart = mart)
av <- results[, "referencedatabase_ensembl", drop=FALSE]
#av <- sapply(av,gsub,pattern="^$",replacement="NA")
#results
#av
#map identifiers
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#iterate vector of ENSEMBL IDs
vec1 <- vector()
vec2 <- vector()
vec3 <- vector()
vec4 <- vector()
for(i in av){
	# query biomart
	results <- getBM(attributes = c("ensembl_gene_id","ensembl_peptide_id","hgnc_symbol","refseq_mrna","entrezgene"), filters = "ensembl_gene_id", values = i, mart = mart)
	vec1 <- c(vec1, results[,2])
	vec2 <- c(vec2, results[,3])
	vec3 <- c(vec3, results[,4])
	vec4 <- c(vec4, results[,5])
}
res1 <- unique(vec1[vec1 != ""])
res2 <- unique(vec2[vec2 != ""])
res3 <- unique(vec3[vec3 != ""])
res4 <- unique(vec4[vec4 != ""])

write.table(res1,paste( args[1], '.ensp', sep=''),row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(res2,paste( args[1], '.hgnc', sep=''),row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(res3,paste( args[1], '.refseq', sep=''),row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(res4,paste( args[1], '.entrez', sep=''),row.names=FALSE,quote=FALSE,col.names=FALSE)

#results[, "symbol"] <- "vec"
#X <- split(results, results$stableidentifier_identifier)
#Y <- lapply(seq_along(X), function(results) as.data.frame(X[[results]])[,5:5])




