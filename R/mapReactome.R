# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#TODO what arguments? What is the true input? How compicated should be filtering?
mapReactome <- function() {
    (nmLipidomics <- read.table("C:/HOME/ONIONpackage/ONION/R/nm-lipidomics.txt",header=TRUE))
    (chebi <- nmLipidomics[1])
    (chebiFiltered <- chebi["ChEBI"])
    (chebiFiltered <- chebi[chebi$ChEBI > 20000 & chebi$ChEBI < 30000,])
    chebiFiltered <- as.data.frame(chebiFiltered)
    colnames(chebiFiltered) <- c("ChEBI")
    chebiFiltered
}

#TODO Reactome API.
mapReactome()


#Original pseudo code API:
#
#data.frame clusterSmallMolecules( file ) {
#    //file - nm-lipidomics.txt (tylko male czasteczki)
#    //K klastrów (Rodzice Chebi)
#    //  H1    H2
#    //1 123   ABC
#    //2 231   DEF
#    //NIE do pliku, jako obiekt.
#    //output - fa_mapping.txt (dzieci, rodzice) w ID z Chebi i chebi parent by ontology
#}

source("R/reactomeAPI.R")
source("R/chebiAPI.R")

clusterSmallMolecules <- function(pathToFile, header=TRUE) {
    smallMolecules <- read.table(pathToFile,header)
    smallMoleculesParents <- lapply(as.list(smallMolecules$ChEBI), function(listElement) {
        #strsplit(listElement$chebiId[1], ":")[[1]][2]

        if (checkIfPathwayIdExistsForChEBIId(listElement)) {
            listElement
        } else {
            listOfChEBIIdsOfOntologyParents <- getListOfChEBIIdsOfOntologyParents(listElement)
            bestParentId <- NA
            for (i in 1:length(listOfChEBIIdsOfOntologyParents)) {
                j <- listOfChEBIIdsOfOntologyParents[[i]]
                if (checkIfPathwayIdExistsForChEBIId(j)) {
                    bestParentId <- j
                    break
                }
            }
            bestParentId
        }
    })
    chEBIIdsToChEBIParentsIDs <- data.frame(ChEBI=as.character(as.list(smallMolecules$ChEBI)), bestChEBIParents=as.character(smallMoleculesParents))
    chEBIIdsToChEBIParentsIDs
}

#Go to chebi, and find at least good identifier in ontology.
#Go to http://www.reactome.org/download/current/ChEBI2Reactome.txt and get Pathways Ids related to CHebi.
#http://www.reactome.org/cgi-bin/instancebrowser?DB=gk_current&ID=6816043&  -  here you can obtain Genes, Chebis and etc.
