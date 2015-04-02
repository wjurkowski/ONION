options(digits = 3, quote = F, scientific = F, stringsAsFactors = F)
library(igraph)
args <- commandArgs(trailingOnly = TRUE)

#read data
t = read.table(args[1], header = F, sep = "\t")
g = graph.edgelist(as.matrix(t[,1:2]))

#betweenness
bt=betweenness(g)
lnBtw=(log2(bt))
lnBtw[lnBtw=="-Inf"] <- 0
n=V(g)$name
#degree
degi=degree(g, mode=c("in"))
dego=degree(g, mode=c("out"))
#closeness
clo=closeness(g)
#eccentricity
ecc=degree(g, mode=c("all"))

#print results 
results<-data.frame(n,degi,dego,lnBtw,ecc,clo)
n_n<-"Node"
degi_n<-"Degree_In"
dego_n<-"Degree_Out"
lnBtw_n<-"ln(BC)"
ecc_n<-"Eccentricity"
clo_n<-"Closeness"
names(results) <-c(n_n,degi_n,dego_n,lnBtw_n,ecc_n,clo_n)
write.table(results, file="NetworkProperties.txt", row.names=F, col.names=T, quote=F, sep="\t")

#global properties i Rmd file

#pa <- get.shortest.paths(g, 5, 9)[[1]]
