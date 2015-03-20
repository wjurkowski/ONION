###################################################
##find interactions for small number of proteins/genes
#define IDs
options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)
v <- read.delim(args[1], header=FALSE)
sp<-as.numeric(args[2])

library(STRINGdb)
string_db <- STRINGdb$new( version="9_1", species=sp,
                           score_threshold=0, input_directory="" )

#first neighbors
first <- vector()
for(i in v){
  a = string_db$mp(i)
  #get neighbourhs
  wyn<-string_db$get_neighbors(a)
  first<-c(first,wyn)
}
uf<-unique(first)

second <- vector()
for(k in uf){
  b = string_db$mp(k)
  #get neighbourhs
  wyn<-string_db$get_neighbors(b)
  second<-c(second,wyn)
}
us<-unique(second)
results<-unique(c(uf,us))

write.table(results,"extended.string",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

#get interactions between listed genes
#string_db$get_interactions( c(tp53, atm) )
#get pubmed Ids of papers supporting interaction between 
#string_db$get_pubmed_interaction( tp53, atm )
