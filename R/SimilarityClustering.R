#library(ChemmineR)
#library(ChemmineOB)


# this function contains preparations for clustering and they are: 1. checking for duplicates (unique_ids),
# 2. printing molecules' IDs as they are named in sdf file, 3. Creating data frame (AllIds) consisting of molecules'
# IDs (these mentioned earlier), it is data frame with one column called ID and number of rows equal to number of molecules,
# 4. Creating counter that is equal to number of rows (number of molecules) in AllIds data frame, counter will be used later in
# while loop, 5. Printing counter (just for checking that everything is ok, that every molecule is in AllIds data frame),
# 6. Creating fingerprints for every molecule (fpset), 7. Print fpset (againf for checking that every molecule has fingerprint),
# 8. Creating empty list (cluster) where created clusters will be stored, each cluster as different component, for
# example cluster[1] will give us names of molecules in cluster number one, custer[2] in cluster number 2 etc.

#TODO: Refactor '<<-', global variable goes to env. Could be some problems with namespaces.
#loading sdf file (containg many sdf molecules).
#SDFFileURI - web address to sdf file/local sdf file.
prep <- function (SDFFileURI) {
    sdfset <- read.SDFset(SDFFileURI)
    unique_ids <- makeUnique(sdfid(sdfset))
    print("IDs are")
    print(cid(sdfset))
    AllIds <<- data.frame(ID=cid(sdfset))
    counter <<- nrow(AllIds)
    cat(paste("number of molecules", counter))
    cat("\n")
    fpset <<- fingerprintOB(sdfset,"FP2")
    print(fpset)
    cluster <<- list()
}

#TODO: Add to public API.
#Check it in tests, here it fails build.
#prep()


# clus () is function for clustering, while loop is used, if counter which is number of rows of AllIds data frame is greater
# than zero then loop is executed. Steps are: 1. Actual clustering using fingerprints created before, method that is used is
# based on Tanimoto index, cutoff can be adjusted to the needs. This function takes the first molecule (from the first row)
# in AllIds data frame (using its fingerprint) and compare it to the rest of molecules that are still in AllIds data frame.
# Results (names of molecules which form the cluster) are stored in 'a' object. 2. Creating data frame called 'b' from
# object 'a'. It consist of one column called ID and as as much rows as there are molecules clustered this time.
# 3. Next step is comparing two data frames: 'AllIds' and 'b' (with clustered molecules) and checking which molecules are present
# in both one data frames - these have FALSE flag and molecules present in only one data frame have TRUE flag.
# These results are stored in 'del' object. 4. Next step is updating of AllIds data frame by deleting these molecules that
# have FALSE flag in 'del' object (becuse they are alraedy clustered). Now in 'AllIds' data frame will be only these molecules
# which are still not clustered. 5. Updating of 'counter' because now 'AllIds' has less rows than before. When 'counter' reaches
# zero (when there will be no molecules in AllIds) while loop won't be excetuted. 6. Updating of 'cluster' which is a list.
# Each new cluster will be added as a  next, new component of this list. Length of this list will be equal to the numbers of
# created clusters.

#TODO: Refactor '<<-', global variable goes to env. Could be some problems with namespaces.
clus <- function() {
    while(counter > 0 ) {
        a <- rownames(as.matrix(fpSim(fpset[as.character(AllIds[1,])], fpset, method="Tanimoto", cutoff=0.3)))
        b <<- data.frame(ID=a)
        del <- (!c(AllIds$ID %in% b$ID))
        AllIds <- AllIds[del,]
        AllIds <- data.frame(ID= AllIds)
        counter <<- nrow(AllIds)
        cluster[[length(cluster)+1]] <<- list(b)
    }
}

#TODO: Add to public API.
#Check it in tests, here it fails build.
#clus()

#checking number of clusters created
amount <- function () {
    print("number of clusters")
    length(cluster)
}

#Add to public API.
#Check it in tests, here it fails build.
#amount()

#This method will cluster every molecule
#,but often one molecule will be found in several clusters.











