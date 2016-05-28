
# NEW PUBLIC API:
clasterUsingOntology <- function(pathToFile, header=TRUE, ontologyRepresentatnion) {
    baseData <- read.table(pathToFile, header)
    ontologyDataFrame <- ontologyRepresentatnion(baseData)
    ontologyDataFrame
}


# NEW PUBLIC API:
mapReactomePathwaysUnderOrganism <- function(chebiOntologyIds, organismTaxonomyId='9606') {
    chebiIdsToEnsembleIds <- ldply(.data = chebiOntologyIds$ontologyId, .fun = function(vectorElement) {
        pathwayIds <- ReactomeAPI::getPathwaysIdsForChebiUnderOrganism(vectorElement, taxonIdToReactomeCodes[[organismTaxonomyId]]$speciesCode)
        ensembleIds <- ReactomeAPI::getEnsemblIdsForPathwayIds(pathwayIds)
        gensSymbols <- getSymbolsBaseOnEnsemblGensIdsUsingMyGenePackage(ensembleIds, organismTaxonomyId = organismTaxonomyId)
        chebiIdToEnsembleIds <- data.frame('chebiId' = as.character(vectorElement)
                                           , 'ensembleIds' = I(list(ensembleIds))
                                           , 'reactomeIds' = I(list(pathwayIds))
                                           , 'gensSymbols' = I(list(gensSymbols)))
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
        stringGensSymbols <- character(length = 0)
        if (0 == length(dfElement$ensembleIds[[1]])) {
        } else {
            proteinIds <- getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(
                dfElement$ensembleIds, organismTaxonomyId = stringOrganismId
            )
            stringId1 <- string_db$mp(proteinIds)
            returnNeighbourVector <- string_db$get_neighbors(stringId1)
            ensembleIdsFromStringDb <- mapFromStringIdsToEnsembleIds(returnNeighbourVector)
            stringGensSymbols <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDb, organismTaxonomyId = stringOrganismId)
        }
        dffff <- data.frame('chebiId' = dfElement$chebiId, 'ensembleIds' = dfElement$ensembleIds[1],
                            'stringIds' = I(list(unique(returnNeighbourVector)))
                            , 'stringGensSymbols' = I(list(unique(stringGensSymbols))) )
        dffff
    })
    dfWithString
}


# NEW API
mapFromStringIdsToEnsembleIds <- function(vactofOfStringIds) {
    ensembleIds <- laply(vactofOfStringIds, .fun = function(vectorElement){
        ensemblePeptideId <- strsplit(vectorElement, "[.]")[[1]][2]
        ensemblePeptideId
    })
    ensembleIds
}


# NEW API.
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


# NEW API.
getSymbolsBaseOnEnsemblGensIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    # genes$symbol
    # genes$ensembl.protein
    additionalInformationBaseOnEnsemblGenId <- queryMany(gensIdsVector, fields = c("symbol","ensembl.protein"),
                                                         species = organismTaxonomyId)
    equivalentEnsemlProteinsIdsVector <- unlist(
        additionalInformationBaseOnEnsemblGenId$symbol[!is.na(additionalInformationBaseOnEnsemblGenId$symbol)]
    )
    equivalentEnsemlProteinsIdsVector <- as.character(equivalentEnsemlProteinsIdsVector)
    equivalentEnsemlProteinsIdsVector
}


# NEW API.
getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    additionalInformationBaseOnEnsemblPeptidId <- queryMany(
        gensIdsVector, scopes = 'ensemblprotein'
        , fields = c("symbol","ensembl.protein"), species = organismTaxonomyId
    )
    equivalentEnsemlProteinsIdsVector <- unlist(
        additionalInformationBaseOnEnsemblPeptidId$symbol[!is.na(additionalInformationBaseOnEnsemblPeptidId$symbol)]
    )
    equivalentEnsemlProteinsIdsVector <- as.character(equivalentEnsemlProteinsIdsVector)
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


# NEW API CCA
makeCanonicalCorrelationAnalysis <- function(xNamesVector, yNamesVector, pathToFileWithXData, pathToFileWithYData) {

    # Where XData = transcriptomicsData and YData = lipidomicsData.
    XData <- readWithoutDuplicates(pathToFileWithXData)
    YData <- readWithoutDuplicates(pathToFileWithYData)

    transposedXData <- as.data.frame(t(XData))
    transposedYData <- as.data.frame(t(YData))

    interX <- intersect(colnames(transposedXData), xNamesVector)
    interY <- intersect(colnames(transposedYData), yNamesVector)

    X <- transposedXData[as.character(interX)]
    Y <- transposedYData[as.character(interY)]

    cca.fit <- NULL

    if (!length(X) || !length(Y)) {
        print("CCA is not possible.")
    } else {
        cca.fit <- yacca::cca(X, Y)
    }
    cca.fit
}

# NEW PUBLIC API
plotCanonicalCorrelationAnalysisResults <- function(ccaResults, x.name = "xLabel", y.name = "yLabel") {
    helio.plot(ccaResults, x.name = "xLabel", y.name = "yLabel")
}


# TODO: Refactor PLS.
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