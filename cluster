library(ChemmineR)
library(ChemmineOB)

#loading sdf file (containg many sdf molecules)
sdfset <- read.SDFset("web address to sdf file/local sdf file")

prep <- function () {
   unique_ids <- makeUnique(sdfid(sdfset))
   print("IDs are")
   print(cid(sdfset))
   AllIds <<- data.frame(ID=cid(sdfset))
   counter <<- nrow(AllIds)
   cat(paste("number of molecules", counter))
   cat("\n")
   fpset <<- fingerprintOB(sdfset,"FP2")
   print(fpset)
   k <<- 1
   cluster <<- list()
}

prep()


clus <- function() {
  while(counter > 0 ) {
    a <- rownames(as.matrix(fpSim(fpset[as.character(AllIds[1,])], fpset, method="Tanimoto", cutoff=0.3)))
    b <<- data.frame(ID=a)
    del <- (!c(AllIds$ID %in% b$ID))
    AllIds <- AllIds[del,]
    AllIds <- data.frame(ID= AllIds)
    counter <<- nrow(AllIds)
    cluster[[length(cluster)+1]] <<- list(b)
    if (counter == 0) break
    k <<- k+1
  }
}

 clus()
 
 # checking number of clusters created
amount <- function () {
  print("number of clusters")
  print(k)
}

amount()
 
 # This method will cluster every molecule but often one molecule will be found in several clusters. 











