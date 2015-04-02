options(digits = 3, quote = F, scientific = F, stringsAsFactors = F)
library(igraph)
args <- commandArgs(trailingOnly = TRUE)

#copy content from script developed in nomo





#read data
t = read.table(args[1], header = F, sep = "\t")

#clusterProfiler
