#R

# read the hypergeom function
source("Enrichment_test.r")

# function that extracts the significant results only
# Coff: cut off value (0.05)
Summary=function(result, Coff=0.05) {
  Sum=result[result$FDR<Coff,][,c(1:7, 10)]
  Sum=Sum[order(Sum$P),]
  return(Sum)
}

# function that counts the nuber of DB IDs in that the dataset was enriched in
# result: result data.frame of the HyperGeomFDR function (from Enrichment_test.r)
# Coff: cut off value (0.05)
P_counter=function(result, Coff=0.05) {
  sign=vector()
  sign=dim(result[result$P<Coff,])[1]
  sign=c(sign, dim(result[result$P_adj_Bonf<Coff,])[1])
  sign=c(sign, dim(result[result$P_adj_BH<Coff,])[1])
  sign=c(sign, dim(result[result$FDR<Coff,])[1])
  names(sign)=c("P", "P_adj_Bonf", "P_adj_BH", "FDR")
  return(sign)
}

# create a list of DB IDs and associated genes from text file
x=strsplit(readLines("REACTOME_database_dmel.txt"), "[[:space:]]+")
RE_genes=lapply(x, tail, n=-1)
names(RE_genes)=lapply(x, head, n=1)
rm(x)

# background genes
BG_genes=readLines("Genes/kept_genes.txt")
# genes to test
Test_genes1=readLines("Genes/test_genes1.txt")
Test_genes2=readLines("Genes/test_genes2.txt")
Test_genes3=readLines("Genes/test_genes3.txt")

# Hypergeometric test for every DB IDs (categories)
# HyperGeomFDR=function(steps, pool, select, DB, nthreads=4) nthreads=max24
Enr_Test_genes1=HyperGeomFDR(1000, BG_genes, Test_genes1, RE_genes, nthreads=20)
write.table(Summary(Enr_Test_genes1, 0.05), file="Enr_Test_genes1_sum.txt", row.names=F, quote=F, sep="\t")
save.image()
Enr_Test_genes2=HyperGeomFDR(1000, BG_genes, Test_genes2, RE_genes, nthreads=20)
write.table(Summary(Enr_Test_genes2, 0.05), file="Enr_Test_genes2_sum.txt", row.names=F, quote=F, sep="\t")
save.image()
Enr_Test_genes3=HyperGeomFDR(1000, BG_genes, Test_genes3, RE_genes, nthreads=20)
write.table(Summary(Enr_Test_genes3, 0.05), file="Enr_Test_genes3_sum.txt", row.names=F, quote=F, sep="\t")
save.image()

# Number of significant results
No_of_P_0.05s=rbind(P_counter(Enr_Test_genes1), P_counter(Enr_Test_genes2), P_counter(Enr_Test_genes3))
rownames(No_of_P_0.05s)=c("Enr_Test_genes1", "Enr_Test_genes1", "Enr_Test_genes1")
write.table(No_of_P_0.05s, file="No_of_P_0.05s.txt", row.names=T, quote=F, sep="\t")
