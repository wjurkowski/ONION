#TODO:  Try remove source in some way.

# source("R/reactomeAPI.R")

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
# NEW PUBLIC API:
clasterUsingOntology <- function(pathToFile, header=TRUE, ontologyRepresentatnion) {
    baseData <- read.table(pathToFile, header)
    ontologyDataFrame <- ontologyRepresentatnion(baseData)
    ontologyDataFrame
}

#PUBLIC API
clusterSmallMolecules <- function(pathToFile, header=TRUE) {
    smallMolecules <- read.table(pathToFile,header)
    chEBIIdsToParentsAndChildren <- decorateChEBITableByFirstParent(smallMolecules)
    chEBIIdsToParentsAndChildren <- decorateChEBITableByFirstChild(chEBIIdsToParentsAndChildren)
}

#PUBLIC API
#TODO: Think about better aproach to <NA>, NA
OLDmergeChEBIOntologyWithChildFavoring <- function(chEBIOntologyUnion) {
    mergeResuls <- mapply(function(firstChEBIParent, firstChEBIChild) {
        if (firstChEBIChild != "NA") {
            as.character(firstChEBIChild)
        } else if (firstChEBIParent != "NA") {
            as.character(firstChEBIParent)
        } else {
            NA
        }
    }, chEBIOntologyUnion$firstChEBIParent, chEBIOntologyUnion$firstChEBIChild)
    data.frame(ChEBI=chEBIOntologyUnion$ChEBI, firstChildOrParent=mergeResuls)
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
    parentIdsToEnsemblIdsdfPL
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

# NEW PUBLIC API:
mapReactomePathwaysUnderOrganism <- function(chebiOntologyIds, organismTaxonomyId='9606') {
    chebiIdsToEnsembleIds <- ldply(.data = chebiOntologyIds$ontologyId, .fun = function(vectorElement) {
        pathwayIds <- ReactomeAPI::getPathwaysIdsForChebiUnderOrganism(vectorElement, taxonIdToReactomeCodes[[organismTaxonomyId]]$speciesCode)
        chebiIdToEnsembleIds <- data.frame('chebiId' = as.character(vectorElement), 'ensembleIds' = I(list(ensembeIds)))
        chebiIdToEnsembleIds
    })
    chebiIdsToEnsembleIds
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



# NEW PUBLIC API
getStringNeighbours <- function(chebiIdsToReactomePathways, stringOrganismId = 9606, stringDbVersion = "10") {
    string_db <- STRINGdb$new( version = stringDbVersion, species = stringOrganismId)
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[!chebiIdsToReactomePathways$chebiId == 0, ]
    dfWithString <- ddply(.data = chebiIdsToRealReactomePathways, .(chebiId), .fun = function(dfElement) {
        returnNeighbourVector <- character(length = 0)
        if (0 == length(dfElement$ensembleIds[[1]])) {
        } else {
            proteinIds <- ONION::getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(dfElement$ensembleIds,
                                                                                             organismTaxonomyId = stringOrganismId)
            stringId1 <- string_db$mp(proteinIds)
            returnNeighbourVector <- string_db$get_neighbors(stringId1)
        }
        dffff <- data.frame('chebiId' = dfElement$chebiId, 'ensembleIds' = dfElement$ensembleIds[1],
                            'stringIds' = I(list(unique(returnNeighbourVector))))
        dffff
    })
}

#PUBLIC API
#Change DbId. To newest one.
OLDgetStringNeighbours <- function (chEBIIdsToGenesSymbolList, stringOrganismId=9606, stringDbVersion = "10") {
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
        if (checkIfValueCanBeProcessed(listElement)) {
            stringIds <- string_db$mp(listElement)
            stringIds
            stringNeignbours <- string_db$get_neighbors(stringIds)
            stringNeignbours
        } else {
            NA
        }
    })
    dataWithStringNeighbours
}

# TODO: Make not public API.
getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    # genes$symbol
    # genes$ensembl.protein
    additionalInformationBaseOnEnsemblGenId <- queryMany(gensIdsVector, fields = c("symbol","ensembl.protein"),
                                                         species = organismTaxonomyId)
    equivalentEnsemlProteinsIds <- unlist(additionalInformationBaseOnEnsemblGenId$ensembl.protein)
    equivalentEnsemlProteinsIdsVector <- lapply(equivalentEnsemlProteinsIds, function(characterListElement){
        characterListElement
    });
    equivalentEnsemlProteinsIdsVector <- as.character(unlist(equivalentEnsemlProteinsIdsVector))
    equivalentEnsemlProteinsIdsVector
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
constans <- list(
    pseudoClasteringsFilters=list(
        ensemblPeptideId="ensembl_peptide_id",
        ensemblGeneId="ensembl_gene_id"
    )
)

#PUBLIC API
showPseudoClusteringResults <- function(chEBIIdsToListOfStringNeigbours, filter) {
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",  host = "jul2015.archive.ensembl.org")
    chEBIIdsToListOfStringNeigboursInEnsembleIds <- mapFromStringIdsToEquivalentEnsembleId(chEBIIdsToListOfStringNeigbours)
    moreDataAbout <- lapply(chEBIIdsToListOfStringNeigboursInEnsembleIds, function(listElement){
        if (length(listElement) == 0) {
            NA
        } else if (is.null(listElement)) {
            NA
        } else if (is.na(listElement)) {
            NA
        } else {
            additionalInformationBaseOnEnsembleGenId <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "refseq_mrna","hgnc_symbol","entrezgene"), filters = c(filter), values = listElement, mart = mart)
            additionalInformationBaseOnEnsembleGenId
        }
    })
}


mapFromStringIdsToEquivalentEnsembleId <- function(chEBIIdsToListOfStringNeigbours) {
    chEBIIdsToListOfStringNeigboursInEnsembleIds <- lapply(chEBIIdsToListOfStringNeigbours, function(listElement) {
        if (length(listElement) == 0) {
            NA
        } else if (is.null(listElement)) {
            NA
        } else if (is.na(listElement)) {
            NA
        } else {
            vectorOfEnsembleProteinIds <- sapply(listElement, function(vectorElement){
                peptideId <- strsplit(vectorElement, "[.]")[[1]][2]
                if (is.na(peptideId)) {
                    vectorElement
                } else {
                    peptideId
                }
            }, USE.NAMES = FALSE);
            vectorOfEnsembleProteinIds
        }
    })
    chEBIIdsToListOfStringNeigboursInEnsembleIds
}




#PUBLIC API
makeCCAOnData <- function(vectorOfChebiIds, vectorOfHgncSymbols,
                          pathToFileWithTranscriptomicsData, pathToFileWithLipidomicsData) {
    transcriptomicsData <- readWithoutDuplicates(pathToFileWithTranscriptomicsData)
    lipidomicsData <- readWithoutDuplicates(pathToFileWithLipidomicsData)
    X <- as.matrix(t(transcriptomicsData))
    Y <- as.matrix(t(lipidomicsData))
    matchedGensData <- match(vectorOfHgncSymbols, colnames(X))
    matchedChebiData <- match(vectorOfChebiIds, colnames(Y))
    factorOfMatchedGensData <- factor(matchedGensData)
    factorOfMatchedChebiData <- factor(matchedChebiData)
    machedX <- X[,as.numeric(levels(factorOfMatchedGensData))]
    machedY <- Y[,as.numeric(levels(factorOfMatchedChebiData))]

    if (is.numeric(machedX) && is.numeric(machedX)) {
        if (is.matrix(machedX)) {
        } else {
            machedX <- matrix(machedX)
        }
        if (is.matrix(machedY)) {
        } else {
            machedY <- matrix(machedY)
        }
    }

    ccaFromyacca <- tryCatch(
        {
            cca(machedX, machedY)
        },
        error=function(cond) {
            message("ONION - Included CCA can not solve task.")
            message("Original message (yacca):")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning=function(cond) {
            message("ONION - Included CCA present warning.")
            message("Original message (yacca):")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
            message("ONION - CCA (yacca) finished with success.")
        }
    )

    #doubleCCA = list(CCA = ccaFromCCA, yacca = ccaFromyacca)
    #doubleCCA
    ccaFromyacca
}

#PUBLIC API
makeGlobalCCA <- function(clasteredList) {
    result <- lapply(as.list(names(clasteredList)), function(listName) {
        listElement <- clasteredList[[listName]]
        if (length(listElement) == 0) {
            NA
        } else if (is.null(listElement)) {
            NA
        } else if (is.na(listElement)) {
            NA
        } else {
            makeCCAOnData(listName, listElement$hgnc_symbol, "D:/doktorat/repositories/ONION/example/transTest.txt",
                          "D:/doktorat/repositories/ONION/example/nm-lipidomics.txt")
        }
    })
    names(result) <- names(clasteredList)
    result
}

showOnlyValuableResults <- function(globalResults) {
    clearGlobalResult <- list()
    savedNames <- c()
    lapply(names(globalResults), function(index, clearGlobalResult){
        if (is.na(globalResults[[index]])) {
        } else {
            clearGlobalResult <<- append(clearGlobalResult, list(index=globalResults[[index]]))
            savedNames <<- c(savedNames, index)
        }
    }, clearGlobalResult=clearGlobalResult)
    names(clearGlobalResult) <- savedNames
    clearGlobalResult
}

makeYaccaCharts <- function(CCAResults, xLabel="GENES", yLabel="LIPIDS") {
    helio.plot(CCAResults, x.name=xLabel, y.name=yLabel)
}

#Analiza różnicowa, differencial analysis.


getOnlyDataMachedInBothSets <- function(vectorOfMoleculesIds, pathToFileWithExperimentalData) {
    transcriptomicsData <- readWithoutDuplicates(pathToFileWithExperimentalData)
    X <- as.matrix(t(transcriptomicsData))
    matchedGensData <- match(vectorOfMoleculesIds, colnames(X))
    factorOfMatchedGensData <- factor(matchedGensData)
    machedX <- X[,as.numeric(levels(factorOfMatchedGensData))]
    machedXMatrix <- convertNumericToMatrix(machedX)
    machedXMatrix
}

convertNumericToMatrix <- function(matrixOrNumericClass) {
    if (is.matrix(matrixOrNumericClassX)) {
        machedX <- matrixOrNumericClass
    } else {
        machedX <- matrix(matrixOrNumericClass)
    }
    machedX
}

makePLSOnData <- function(vectorOfChebiIds, vectorOfHgncSymbols,
                          pathToFileWithTranscriptomicsData, pathToFileWithLipidomicsData) {
    machedXMatrix <- getOnlyDataMachedInBothSets(vectorOfHgncSymbols, pathToFileWithTranscriptomicsData)
    machedYMatrix <- getOnlyDataMachedInBothSets(vectorOfChebiIds, pathToFileWithLipidomicsData)

    PLSResults <- tryCatch(
        {
            Xmelt <- I(as.matrix(machedXMatrix))
            Ymelt <- I(as.matrix(machedYMatrix))
            realData = list(Xmelt,Ymelt)
            # Coefficiens - what should be, or what type of arguments?
            #, ncomp=as.numeric("10")
            plsr(Y ~ X, data = realData, validation="LOO")
        },
        error=function(cond) {
            message("ONION - PLS fail, it can not solve problem.")
            message("Original message (PLS):")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning=function(cond) {
            message("ONION - PLS fail, it can not solve problem.")
            message("Original message (PLS):")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
            message("ONION - PLS finished with success.")
        }
    )

    PLSResults
}

makePLSCharts <- function(PLS) {
    png("PLS_loadings.png", width = 640, height = 480)
        par(mfrow = c(2,2))
        biplot(PLS, which = "x") # Default
        biplot(PLS, which = "y")
        biplot(PLS, which = "scores")
        biplot(PLS, which = "loadings")
    dev.off()
}