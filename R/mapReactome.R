source("R/reactomeAPI.R")
source("R/chebiAPI.R")
library(biomaRt)
library(STRINGdb)

#Original pseudo code API (OPCAPI):
#
#data.frame clusterSmallMolecules( file ) {
#    //file - nm-lipidomics.txt (tylko male czasteczki)
#    //K klastr?w (Rodzice Chebi)
#    //  H1    H2
#    //1 123   ABC
#    //2 231   DEF
#    //NIE do pliku, jako obiekt.
#    //output - fa_mapping.txt (dzieci, rodzice) w ID z Chebi i chebi parent by ontology
#}

#PUBLIC API
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
            bestChildId <- NA
            for (i in 1:length(listOfChEBIIdsOfOntologyChildren)) {
                j <- listOfChEBIIdsOfOntologyChildren[[i]]
                if (checkIfPathwayIdExistsForChEBIId(j)) {
                    bestChildId <- j
                    break
                }
            }
            bestChildId
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

#PUBLIC API
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

#PUCLIC API
#TODO: Add specification for organism.
mapReactomePathways <- function(chEBIToParentsIdDataFrame, organismCode='9606') {
    parentIdsList <- as.list(as.vector(chEBIToParentsIdDataFrame$firstChildOrParent))
    parentIdsToPathways <- lapply(parentIdsList, function(listElement) {
        if (!is.na(listElement)) { #listElement != "NA") {
            listOfUsablePathwaysIds <- getPathwaysIdsForChEBIAndOrganismCode(listElement, taxonIdToReactomeCodes[[organismCode]]$speciesCode)
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

    NoDuplicatesOnVectorElementsAndCleanFlatParentIdToEnsemblIdsList <- removeDuplicatesInListVectors(cleanAndFlatParentIdToEnsemblIdsList)

    ad <- addNames(chEBIToParentsIdDataFrame, NoDuplicatesOnVectorElementsAndCleanFlatParentIdToEnsemblIdsList)
    ad
    #    candidateSet 2975806 //Reactome - Set
    #    entityWithAccessionedSequence 193137 //ENSEMBL, UniProt - Protein
    #    simpleEntity 113533 //ChEBI - Chemical Compound
    #    complex 5580259 //Reactome - Complex
}

#Static hashmap taxonId <-> (reactome$speciesName, reactome$spaciesCode) Code is used in resouce files.
taxonIdToReactomeCodes <- new.env()
taxonIdToReactomeCodes[['9606']] <- list(speciesName='Homo sapiens', speciesCode='HSA')
taxonIdToReactomeCodes[['3702']] <- list(speciesName='Arabidopsis thaliana', speciesCode='ATH')
taxonIdToReactomeCodes[['9913']] <- list(speciesName='Bos taurus', speciesCode='BTA')
taxonIdToReactomeCodes[['6239']] <- list(speciesName='Caenorhabditis elegans', speciesCode='CEL')
taxonIdToReactomeCodes[['9615']] <- list(speciesName='Canis familiaris', speciesCode='CFA')
taxonIdToReactomeCodes[['7955']] <- list(speciesName='Danio rerio', speciesCode='DRE')
taxonIdToReactomeCodes[['44689']] <- list(speciesName='Dictyostelium discoideum', speciesCode='DDI')
taxonIdToReactomeCodes[['7227']] <- list(speciesName='Drosophila melanogaster', speciesCode='DME')
taxonIdToReactomeCodes[['9031']] <- list(speciesName='Gallus gallus', speciesCode='GGA')
taxonIdToReactomeCodes[['10090']] <- list(speciesName='Mus musculus', speciesCode='MMU')
taxonIdToReactomeCodes[['1773']] <- list(speciesName='Mycobacterium tuberculosis', speciesCode='MTU')
taxonIdToReactomeCodes[['4530']] <- list(speciesName='Oryza sativa', speciesCode='OSA')
taxonIdToReactomeCodes[['5833']] <- list(speciesName='Plasmodium falciparum', speciesCode='PFA')
taxonIdToReactomeCodes[['10116']] <- list(speciesName='Rattus norvegicus', speciesCode='RNO')
taxonIdToReactomeCodes[['4932']] <- list(speciesName='Saccharomyces cerevisiae', speciesCode='SCE')
taxonIdToReactomeCodes[['4896']] <- list(speciesName='Schizosaccharomyces pombe', speciesCode='SPO')
taxonIdToReactomeCodes[['8364']] <- list(speciesName='Xenopus tropicalis', speciesCode='XTR')
taxonIdToReactomeCodes[['59729']] <- list(speciesName='Taeniopygia guttata', speciesCode='TGU')
taxonIdToReactomeCodes[['9823']] <- list(speciesName='Sus scrofa', speciesCode='SSC')

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

removeDuplicatesInListVectors <- function(listOfVectors) {
    listOfVectorsWithoutDuplicates <- lapply(listOfVectors, function(listElement){
        unique(x = listElement)
    })
    listOfVectorsWithoutDuplicates
}

addNames <- function(ChEBIList, flattenedList) {
    for (i in 1:length(flattenedList)) {
        names(flattenedList)[i] <- ChEBIList[i,]$ChEBI
    }
    flattenedList
}

removeEmptyElementsOnListOfCharsVectors <- function(listOfVectors) {
    listOfVectorsWithoutEmptyChars <- lapply(listOfVectors, function(listElement){
        clearEnsemblePeptideIds <- listElement[listElement != ""]
        clearEnsemblePeptideIds
    })
    listOfVectorsWithoutEmptyChars
}

#PUBLIC API
#Change DbId. To newest one.
getStringNeighbours <- function (chEBIIdsToGenesSymbolList, stringOrganismId=9606, stringDbVersion = "10") {
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",  host = "jul2015.archive.ensembl.org")
    moreDataAbout <- lapply(chEBIIdsToGenesSymbolList, function(listElement){
        if (is.null(listElement)) {
            "NA"
        } else {
            additionalInformationBaseOnEnsembleGenId <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "refseq_mrna","hgnc_symbol","entrezgene"), filters = c("ensembl_gene_id"), values = listElement, mart = mart)
            #g[['28125']][g[['28125']] != ""]
            toclean <- additionalInformationBaseOnEnsembleGenId$ensembl_peptide_id
            toclean
            clearEnsemblePeptideIds <- toclean[toclean != ""]
            clearEnsemblePeptideIds
        }
    })
    moreDataAboutWithoutEmptyElements <- removeEmptyElementsOnListOfCharsVectors(moreDataAbout)
    noDuplicates <- removeDuplicatesInListVectors(moreDataAboutWithoutEmptyElements)
    noDuplicates
    string_db <- STRINGdb$new( version=stringDbVersion, species=stringOrganismId, score_threshold=0, input_directory="")
    dataWithStringNeighbours <- lapply(noDuplicates, function(listElement){
        stringIds <- string_db$mp(listElement)
        stringIds
        stringNeignbours <- string_db$get_neighbors(stringIds)
        stringNeignbours
    })
    dataWithStringNeighbours
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


#PUBLIC API
showPseudoClusteringResults <- function(chEBIIdsToListOfStringNeigbours) {
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",  host = "jul2015.archive.ensembl.org")
    chEBIIdsToListOfStringNeigboursInEnsembleIds <- mapFromStringIdsToEnsembleProteinId(chEBIIdsToListOfStringNeigbours)
    moreDataAbout <- lapply(chEBIIdsToListOfStringNeigboursInEnsembleIds, function(listElement){
        if (is.null(listElement)) {
            "NA"
        } else if (length(listElement) == 0 ) {
            "NA"
        } else {
            additionalInformationBaseOnEnsembleGenId <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "refseq_mrna","hgnc_symbol","entrezgene"), filters = c("ensembl_peptide_id"), values = listElement, mart = mart)
            additionalInformationBaseOnEnsembleGenId
        }
    })
}


mapFromStringIdsToEnsembleProteinId <- function(chEBIIdsToListOfStringNeigbours) {
    chEBIIdsToListOfStringNeigboursInEnsembleIds <- lapply(chEBIIdsToListOfStringNeigbours, function(listElement) {
        vectorOfEnsembleProteinIds <- sapply(listElement, function(vectorElement){
            peptideId <- strsplit(vectorElement, "[.]")[[1]][2]
            peptideId
        }, USE.NAMES = FALSE)
        vectorOfEnsembleProteinIds
    })
    chEBIIdsToListOfStringNeigboursInEnsembleIds
}
