library(ChemmineR)
library(ChemmineOB)

# loading file ( read.SDFset imports one or many molecules from an SD/MOL file and stores it in an
# SDFset container. Supports both the V2000 and V3000 formats.)
sdfset <- read.SDFset(path to the local SD file)

# this function is just for printing molecules' IDs and names and is optional (not necessary to run).
# This function takes as an argument sdfset created/loaded earlier. 
showNames <- function (sdfset) {
  cat("IDs are: \n") 
  print(cid(sdfset))
  cat("\n")
  cat("molecules' names are: \n") 
  print(sdfid(sdfset)) 
}

showNames(sdfset)


# This function contains preparations for clustering and they are: 
# 1. Checking for duplicates (unique_ids), 
# 2. Identification of invalid SDFs in SDFset objects and removing invalid SDFs, if there are any, 
# 3. Creating two identical data frames (AllIds and AllIds2) - one for each of the following functions (clusterWithCutoff and 
# makeSimilarityMatrix) - consisting of molecules' IDs (as mentioned earlier), they are data frames with 
# only one column called ID and number of rows equal to number of molecules,
# 4. Creating two identical counters (counter and counter2) for the same reasons as two AllIds data frames. 
# These counters are equal to number of rows (number of molecules) in AllIds/AllIds2 data frames, 
# counters will be used later in both while loops, 
# 5. Printing counter (just for checking that everything is ok and the number is all right), printig just one
# counter as the second is identical
# 6. Creating fingerprints for every molecule (fpset)
# 7. Print fpset (againf for checking that every molecule has fingerprint).
# This function takes as an argument sdfset created/loaded earlier. 


prep <- function (sdfset) {
  print(sdfset) # to check if loaded sdf file contains as much molecules as it should be
  unique_ids <- makeUnique(sdfid(sdfset))
  valid <- validSDF(sdfset) #  
  sdfset <- sdfset[valid] # 
  AllIds <<- data.frame(ID=cid(sdfset))
  AllIds2 <<- data.frame(ID=cid(sdfset))
  counter <<- nrow(AllIds)
  counter2 <<- nrow(AllIds2)
  print(counter)
  cat(paste("number of molecules:", counter, "\n")) 
  fpset <<- fingerprintOB(sdfset,"FP2")
  print(fpset)
}

prep(sdfset)

# setting a cutoff for clustering (by default it's 0.3 but can be adjusted to the needs)
cutoff <- 0.3


# clusterWithCutoff() is a function for similarity clustering with a cutoff. It takes previously created fingerprints for 
# every molecule and it does the clustering by using Tanimoto method on these fingerprints with the cutoff set on 0.3 (on default) which
# means that every molecule that after clustering mentioned above has score under 0.3 won't be a part of this cluster. 
# This function keeps on repeating until there is no more molecules that are not clustered (AllIds is empty, counter = 0).
# Below are the steps which make up this function and they are as follows:
# 1. Creating empty list (cluster) where created clusters will be stored, each cluster as different component of this list,
# for example cluster[1] will give us IDs of molecules in cluster number one, custer[2] in cluster number 2 etc.
# Calling "cluster" after running this function will provide user with list of created clusters.
# 2. While loop is used with value of the counter as condtion.  If the value of the counter (which is number of rows of 
# AllIds data frame) is greater than zero then loop is executed. It keeps on repeating until counter = 0
# which means there are no more molecules to be clustered and loop is stopped. 
# 3. Inside while loop takes place an actual clustering. As mentioned before clustering is done by using method 
# based on Tanimoto index on previously created fingerprints with a set cutoff (by default it's 0.3 ).
# It takes the first molecule (from the first row) from AllIds data frame, using its fingerprint and compare it 
# to the rest of molecules' fingerprints that are still in AllIds data frame. Results of that action 
#  are stored in 'a' object as IDs of clustered molecules.  
# 4. Creating data frame called 'b' from object 'a'.  This new data frame consist of one column called ID and 
# as as much rows as there are molecules clustered this time.
# 5. Next step is comparing two data frames: 'AllIds' and 'b' (with clustered molecules) and checking which molecules are 
# present in both data frames - these will get FALSE flag and molecules present in only one of these data frames will get 
# TRUE flag. These results (FALSEs and TRUEs) are stored in 'del' object. 
# 6. Next step is updating of AllIds data frame by deleting these molecules that have FALSE flag in 'del' object 
# (becuse they are already clustered). Now in 'AllIds' data frame will be only these molecules
# which are still not clustered. 
# 7.Updating of 'counter' because now 'AllIds' has less rows than before. When 'counter' reaches
# zero (when there will be no molecules in AllIds left) while loop won't be excetuted.
# 8. Updating of 'cluster' which is a list. Each new cluster will be added as a next, new component of 
# this list. Length of this list will be equal to the numbers of created clusters.
# This function takes as an argument fpset (object with molecules' fingerprints) created earlier.
# This function will cluster all molecules but one molecule can be found in only one cluster because 
# if the molecule is already clustered then it's removed from pool of molecules still waiting to be clustered.

clusterWithCutoff <- function(fpset) {
  cluster <<- list()
  while(counter > 0 ) {
    a <- rownames(as.matrix(fpSim(fpset[as.character(AllIds[1,])], fpset, method="Tanimoto", cutoff= cutoff)))
    b <- data.frame(ID=a)
    del <- (!c(AllIds$ID %in% b$ID))
    AllIds <- AllIds[del,]
    AllIds <- data.frame(ID= AllIds)
    counter <- nrow(AllIds)
    cluster[[length(cluster)+1]] <<- list(b)
  }
  return(cluster)
}

clusterWithCutoff(fpset)


# This function is just for checking the number of created clusters. It's not necessary to run. 
# This function takes as an argument object "cluster" created earlier.
amount <- function(cluster) {
  print("number of clusters")
  length(cluster)
}

amount(cluster)


# makeSimilarityMatrix() is a function that creates a matrix with similarity values between any two molecules stored
# in AllIds2 data frame (which is created from "sdfset" object and is identical to Allids data frame). 
# Like clusterWithCutoff function this one does the clustering by using Tanimoto method on molecules' fingerprints 
# but without cutoff (cutoff is set to 0.0). Thanks to that as a result we get similarity matrix with 
# n-rows and n-columns (square) where n is a number of molecules, so there isn't a pair of molecules that don't
# have similarity value calculated. Similarity between the two identical molecules are equal to 1 and makes
# the diagonal of this matrix. 
# This function keeps on repeating until there is no more molecules that are not have their similarity
# with other molecules computed (AllIds is empty, counter = 0).
# Below are the steps which make up this function and they are as follows:
# 1. Creating empty matrix (values) where created data frames will be stored, one data frame per molecule,
# data frame fo every molecule. The numbers of rows and columns are the same as number of molecules.
# Each new data frame with the computed values between given molecule and the rest of molecules (and between
# molecule identical to the one given too) will become a new column in created matrix. So every new data frame
# will be added to the mtarix as new column.
# 2. While loop is used with value of the counter2 as condtion. If the value of the counter2 (which is number of rows of 
# AllIds2 data frame) is greater than zero then loop is executed. It keeps on repeating until counter = 0
# which means there are no more molecules which have not been compared with the others. 
# 3. Inside while loop takes place an actual silimarity computing and matrix creating. As mentioned before clustering is 
# done by using method based on Tanimoto index on previously created fingerprints without a cutoff (cutoff = 0.0).
# This time cutoff  must be equal to 0 and can't be changed because otherwise we won't get full similarity
# matrix (between any two molecules).
# It takes the first molecule (from the first row) from AllIds2 data frame, using its fingerprint and compare it 
# to the rest of molecules' fingerprints and to itself too. Results of that action which are similarity scores for 
# any pair of molecules  are stored in 'd' object. 
# 4. Setting molecules' IDs as columns names. 
# 4. Updating AllIds2 data frame by changing the first molecule, so the next molecule will be new number 
# one in this data frame. 
# 5.Updating of 'counter2' because now 'AllIds2' has less rows than before. When 'counter2' reaches
# zero (when there will be no molecules in AllIds2 left) while loop won't be excetuted.
# 6. Creating a new object "n" wich is equal to the number of rows in created matrix. It is necessary to do
# because it will be used as condition in later if/else loop.
# 7. If/else loop where if n == 1 then one new column will be added to the created, still empty matrix. 
# It will be done only once at the beginning because after adding new column "n" won't be equal to 1 any longer.
# Then else will be executed where function "cbind" is used - it binds new data frame as new column to the
# already existing matrix. 
# This function takes as an argument fpset (object with molecules' fingerprints) created earlier.


makeSimilarityMatrix <- function (fpset) {
values <<- matrix()

  while (counter2 > 0) {
    d <- as.data.frame(fpSim(fpset[as.character(AllIds2[1,])], fpset, sorted=FALSE, method="Tanimoto", cutoff=0.0))
    colnames(d) <- as.character(AllIds2[1,])
    AllIds2 <- data.frame(AllIds2[-1, ])
    counter2 <- nrow(AllIds2)
    n <- nrow(values)
    if (r == 1){
      values <<- d
    }
    else {
      values <<- cbind(values,d)
    }
  }
  return(values)
}

makeSimilarityMatrix(fpset)

# full, created matrix can be exported to the .csv file
write.csv(values, "SimilarityMatrix.csv")









