source("R/reactomeAPI.R")
source("R/chebiAPI.R")

#Original pseudo code API (OPCAPI):
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

clusterSmallMolecules <- function(pathToFile, header=TRUE) {
    smallMolecules <- read.table(pathToFile,header)
    chEBIIdsToParentsAndChildren <- decorateChEBITableByFirstParent(smallMolecules)
    chEBIIdsToParentsAndChildren <- decorateChEBITableByFirstChild(chEBIIdsToParentsAndChildren)
}

decorateChEBITableByFirstChild <- function(ChEBITable) {
    smallMoleculesChildren <- lapply(as.list(as.character(ChEBITable$ChEBI)), function(listElement) {
        if (checkIfPathwayIdExistsForChEBIId(listElement)) {
            listElement
        } else {
            listOfChEBIIdsOfOntologyChildren <- getListOfChEBIIdsOfOntologyChildren(listElement)
            bestParentId <- NA
            for (i in 1:length(listOfChEBIIdsOfOntologyChildren)) {
                j <- listOfChEBIIdsOfOntologyChildren[[i]]
                if (checkIfPathwayIdExistsForChEBIId(j)) {
                    bestParentId <- j
                    break
                }
            }
            bestParentId
        }
    })
    ChEBITableAndFirstChEBIChildUnion <- data.frame(ChEBITable, firstChEBIChild=as.character(smallMoleculesChildren))
}

decorateChEBITableByFirstParent <- function(ChEBITable) {
    smallMoleculesParents <- lapply(as.list(as.character(ChEBITable$ChEBI)), function(listElement) {
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
    ChEBITableAndFirstChEBIParentUnion <- data.frame(ChEBITable, firstChEBIParent=as.character(smallMoleculesParents))
}

#TODO: Think about better aproach to <NA>, NA
margeChEBIOntologyWithChildFavoring <- function(chEBIOntologyUnion) {
    margeResuls <- mapply(function(firstChEBIParent, firstChEBIChild) {
        if (firstChEBIChild != "NA") {
            as.character(firstChEBIChild)
        } else if (firstChEBIParent != "NA") {
            as.character(firstChEBIParent)
        } else {
            NA
        }
    }, chEBIOntologyUnion$firstChEBIParent, chEBIOntologyUnion$firstChEBIChild)
    data.frame(ChEBI=chEBIOntologyUnion$ChEBI, firstChildOrParent=margeResuls)
}

#Go to chebi, and find at least good identifier in ontology.
#Go to http://www.reactome.org/download/current/ChEBI2Reactome.txt and get Pathways Ids related to CHebi.
#http://www.reactome.org/cgi-bin/instancebrowser?DB=gk_current&ID=6816043&  -  here you can obtain Genes, Chebis and etc.


#OPCAPI:
#
#list[list] mapReactomePathways <- function ( data.frame fa_mapping.txt, ) {
#    //input - fa_mapping.txt
#    //all mapReactome.R file
#    //MY: Find not duplicated PAthways
#    //TODO: find Gens.
#    //smallMolecules[Genes], NOW:smallMolecules[Pathways]
#    //output genesSymbolList smallMolecules[Genes] K[x1...xn]
#}

#TODO: Add specification for organism.
mapReactomePathways <- function(chEBIToParentsIdDataFrame, organismCode) {
    parentIdsList <- as.list(as.vector(chEBIToParentsIdDataFrame$firstChildOrParent))
    parentIdsToPathways <- lapply(parentIdsList, function(listElement) {
        if (!is.na(listElement)) { #listElement != "NA") {
            listOfUsablePathwaysIds <- getPathwaysIdsForChEBIAndOrganismCode(listElement, organismCode)
        } else {
            NA
        }
    })
    parentIdsToPathways
    parentIdsToEnsemblIds <- lapply(parentIdsToPathways, function(listElement){
        abc <- lapply(listElement, function(pathwayId) {
            getEnsemblIdsForPathway(pathwayId)
        })
    })
    parentIdsToEnsemblIds
    cleanParentIdToEnsemblIdsList <- lapply(parentIdsToEnsemblIds, function(listElement) {
        removeEmptyLists(listElement)
    })
    cleanAndFlatParentIdToEnsemblIdsList <- flattenList(cleanParentIdToEnsemblIdsList)

    ad <- addNames(chEBIToParentsIdDataFrame, cleanAndFlatParentIdToEnsemblIdsList)

#    candidateSet 2975806 //Reactome - Set
#    entityWithAccessionedSequence 193137 //ENSEMBL, UniProt - Protein
#    simpleEntity 113533 //ChEBI - Chemical Compound
#    complex 5580259 //Reactome - Complex
}

#TODO move to utils
#priver util methods
removeEmptyLists <- function(dirtyList) {
    recurention <- FALSE
    if (length(dirtyList) == 0) {
    } else {
        for (i in 1:length(dirtyList)) {
            if (length(dirtyList[[i]]) == 0) {
                dirtyList[[i]] <- NULL
                recurention <- TRUE
                break
            }
        }
    }
    if (recurention) {
        dirtyList <- removeEmptyLists(dirtyList)
    }
    dirtyList
}

flattenList <- function(listToFlatten) {
    lapply(listToFlatten, function(listElement){
        if (length(listElement) != 0) {
            listToReturn <- list()
            for (i in 1:length(listElement)) {
                if (length(listElement[[i]]) != 0) {
                    for (j in 1:length(listElement[[i]])) {
                        listToReturn <- c(listToReturn, listElement[[i]][[j]])
                    }
                }
            }
            #as.vector(listToReturn)
            as.character(listToReturn)
            #listToReturn
        }
    })
}

addNames <- function(ChEBIList, flattenedList) {
    for (i in 1:length(flattenedList)) {
        names(flattenedList)[i] <- ChEBIList[i,]$ChEBI
    }
    flattenedList
}

#TODO: NOW
#Use Ensamble2Reactome.txt file
#
#Use this files:
#UniProt to pathways mapping file
#ChEBI to pathways mapping file
#ENSEMBL to pathways mapping file
#
#Translate (Biomart, clusterProfilerBitR, myGene)
#
#
#
#TODO: FUTURE
#http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/pathwayParticipants/109581
#->
#http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/queryById/entityWithAccessionedSequence/197679
#->
#ENSEMBLE,
#
#
#


