
# NEW PUBLIC API:
clusterUsingOntology <- function(chebiIdsDataFrame, rootColumnName, ontologyRepresentatnion) {
    ontologyDataFrame <- ontologyRepresentatnion(baseData = chebiIdsDataFrame, rootColumnName = rootColumnName)
    ontologyDataFrame
}


# NEW PUBLIC API:
mapReactomePathwaysUnderOrganism <- function(chebiOntologyIds, organismTaxonomyId='9606', idsColumnName = 'ontologyId', rootColumnName = 'root') {
    # x <- c("supp", "dose")
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(idsColumnName, rootColumnName)
    }
    chebiIdsToEnsembleIds <- ddply(.data = chebiOntologyIds, columnsUseInIteration, .fun = function(vectorElement) {

        print(as.character(vectorElement[1, c(idsColumnName)]))
        idToCheck <- as.character(strsplit(as.character(vectorElement[1, c(idsColumnName)]), ":")[[1]][2])
        print("idToCheck")
        print(idToCheck)
        if (is.na(idToCheck)) {
            idToCheck <- as.character("0");
            print(idToCheck)
        }
        pathwayIds <- getPathwaysIdsForChebiUnderOrganism(idToCheck, taxonIdToReactomeCodes[[organismTaxonomyId]]$speciesCode)
        ensembleIds <- getEnsemblIdsForPathwayIds(pathwayIds)
        uniProtIds <- getUniProtIdsForPathwayIds(pathwayIds)
        genesSymbolsFromEnsemble <- getSymbolsBaseOnEnsemblGensIdsUsingMyGenePackage(ensembleIds, organismTaxonomyId = organismTaxonomyId)
        genesSymbolsFromUniProt <- getSymbolsBaseOnUniProtIdsUsingMyGenePackage(uniProtIds, organismTaxonomyId = organismTaxonomyId)
        chebiIdToEnsembleIds <- data.frame('ensembleIds' = I(list(ensembleIds)),
                                           'uniProtIds' = I(list(uniProtIds)),
                                           'reactomeIds' = I(list(pathwayIds)),
                                           'genesSymbolsFromEnsemble' = I(list(genesSymbolsFromEnsemble)),
                                           'genesSymbolsFromUniProt' = I(list(genesSymbolsFromUniProt)))
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
getStringNeighbours <- function(chebiIdsToReactomePathways, stringOrganismId = 9606, stringDbVersion = "10", idsColumnName = 'ontologyId', rootColumnName = 'root', listOfEnsembleIdColumnName = 'ensembleIds') {
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(idsColumnName, rootColumnName)
    }
    string_db <- STRINGdb$new( version = stringDbVersion, species = stringOrganismId)
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[!chebiIdsToReactomePathways[idsColumnName] == 0, ]
    dfWithString <- ddply(.data = chebiIdsToRealReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
        returnNeighbourVector <- character(length = 0)
        stringGenesSymbols <- character(length = 0)
        if (0 == length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
        } else {
            proteinIds <- getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(
                dfElement[1,listOfEnsembleIdColumnName], organismTaxonomyId = stringOrganismId
            )
            stringId1 <- string_db$mp(proteinIds)
            returnNeighbourVector <- string_db$get_neighbors(stringId1)
            ensembleIdsFromStringDb <- mapFromStringIdsToEnsembleIds(returnNeighbourVector)
            stringGenesSymbols <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDb, organismTaxonomyId = stringOrganismId)
        }
        dffff <- data.frame('ensembleIds' = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(returnNeighbourVector))),
                            'stringGenesSymbols' = I(list(unique(stringGenesSymbols))) )
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
getSymbolsBaseOnUniProtIdsUsingMyGenePackage <- function(gensIdsVector, organismTaxonomyId) {
    # genes <- getGenes(gensIdsVector, fields = "all")
    # genes$symbol
    # genes$ensembl.protein
    additionalInformationBaseOnEnsemblGenId <- queryMany(gensIdsVector, scopes = 'uniprot', fields = c("symbol"),
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

# NEW PUBLIC API
groupUsingUserDefinition <- function(pathToFileWithGroupDefinition ) {
    maxColLength <- max(count.fields(pathToFileWithGroupDefinition, sep = '\t'))
    model <- read.table(file = pathToFileWithGroupDefinition, header = TRUE, fill = TRUE,
                        stringsAsFactors = FALSE, sep = "\t", strip.white = TRUE)
}



# NEW API CCA
makeCanonicalCorrelationAnalysis <- function(xNamesVector, yNamesVector, XDataFrame, YDataFrame) {

    # Where XData = transcriptomicsData and YData = lipidomicsData.
    XData <- data.frame(XDataFrame[!duplicated(XDataFrame[1]), ], row.names = 1)
    YData <- data.frame(YDataFrame[!duplicated(YDataFrame[1]), ], row.names = 1)

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
        # print("X IS : ")
        # print(X)
        # print("Y IS : ")
        # print(Y)
        # cca.fit <- yacca::cca(X, Y)


        cca.fit <- tryCatch(
            {
                yacca::cca(X, Y)
            },
            error = function(cond) {
                message("ONION - Included CCA can not solve task.")
                message("Original message (yacca):")
                message(cond)
                # Choose a return value in case of error
                return(NA)
            },
            warning = function(cond) {
                message("ONION - Included CCA present warning.")
                message("Original message (yacca):")
                message(cond)
                # Choose a return value in case of warning
                return(NULL)
            },
            finally = {
                message("ONION - CCA (yacca) finished.")
            }
        )
    }
    cca.fit
}


# NEW PUBLIC API
makePermutationTestOnCCA <- function(XDataFrame, YDataFrame, numberOfRowsForTestOnX, numberOfRowsForTestOnY, numberOfIterations = 100, countedCCA) {
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(numberOfRowsForTestOnX)
    print(numberOfRowsForTestOnY)
    # print()
    # print()
    vectorOfXrd <- as.numeric();
    vectorOfYrd <- as.numeric();
    for (i in 1:numberOfIterations) {
        xNV <- as.character(XDataFrame[sample(nrow(XDataFrame), numberOfRowsForTestOnX), ][,1])
        yNV <- as.character(YDataFrame[sample(nrow(YDataFrame), numberOfRowsForTestOnY), ][,1])
        ccaResult <- ONION::makeCanonicalCorrelationAnalysis(xNamesVector = xNV, yNamesVector = yNV, XDataFrame = XDataFrame, YDataFrame = YDataFrame)
        print("++++++++++++++++++++++++++++++++++")
        # print.default(ccaResult)
        if (is.na(ccaResult) || is.null(ccaResult)) {

        } else {
            vectorOfXrd <- c(vectorOfXrd, ccaResult$xrd)
            vectorOfYrd <- c(vectorOfYrd, ccaResult$yrd)
        }
    }

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print.default(countedCCA)
    print("*****************************")
    print(countedCCA$xrd)
    print(countedCCA$yrd)
    print(vectorOfXrd)
    print(vectorOfYrd)
    print(mean(vectorOfXrd))
    print(mean(vectorOfYrd))
    meanOfXrd <- NA;
    meanOfYrd <- NA;
    if (0 != length(vectorOfXrd)) {
        meanOfXrd <- mean(vectorOfXrd);
    }
    if (0 != length(vectorOfYrd)) {
        meanOfYrd <- mean(vectorOfYrd);
    }
    testResult <- list("countedCCAOnX" = countedCCA$xrd,
                       "countedCCAOnY" = countedCCA$yrd,
                       "meanOnX" = meanOfXrd,
                       "meanOnY" = meanOfYrd)
    testResult
}


# NEW PUBLIC API
makeCCAOnGroups <- function(groupsDefinitionDF, mappingDF, leftMappingColumnName = 'root', rightMappingColumnName = 'genesSymbolsFromEnsemble', groupsDataDF, mappingDataDF){
    ddply(.data = groupsDefinitionDF['Molecules'], .(Molecules), .fun = function(dfElement) {
        print("???????????????????")
        print(dfElement)
        rightSideIdsToAnalys <- unlist(strsplit(as.character(dfElement), split = " "));
        print(rightSideIdsToAnalys)
        leftSideIdsToAnalys <- mappingDF[mappingDF[[leftMappingColumnName]] %in% rightSideIdsToAnalys,][[rightMappingColumnName]]
        leftSideIdsToAnalys <- unique(unlist(leftSideIdsToAnalys))

        ccaResults <- ONION::makeCanonicalCorrelationAnalysis(
            xNamesVector = leftSideIdsToAnalys,
            yNamesVector = rightSideIdsToAnalys,
            XDataFrame = mappingDataDF,
            YDataFrame = groupsDataDF)

        print("######################################")
        # print(ccaResults)
        print("---------------------------")
        # print.default(ccaResults)
        #TODO : Use user defined column name instead of symbol. ChEBI column too.
        numberOfRowsForTestOnX <- nrow(mappingDataDF[mappingDataDF$symbol %in% leftSideIdsToAnalys, ])
        numberOfRowsForTestOnY <- nrow(groupsDataDF[groupsDataDF$ChEBI %in% rightSideIdsToAnalys, ])
        parmutationTestResult <- makePermutationTestOnCCA(XDataFrame = mappingDataDF, YDataFrame = groupsDataDF,
                                                          numberOfRowsForTestOnX = numberOfRowsForTestOnX,
                                                          numberOfRowsForTestOnY = numberOfRowsForTestOnY,
                                                          numberOfIterations = 50, countedCCA = ccaResults);

        dfWithCca <- data.frame('right' = I(list(rightSideIdsToAnalys)),
                                'left' = I(list(leftSideIdsToAnalys)),
                                'ccaResults' = I(list(ccaResults)),
                                'ccaPermutationTestResults' = I(list(parmutationTestResult)))
        dfWithCca
    })
}


# NEW PUBLIC API
plotCanonicalCorrelationAnalysisResults <- function(ccaResults, x.name = "xLabel", y.name = "yLabel") {
    helio.plot(ccaResults, x.name = "xLabel", y.name = "yLabel")
}


# NEW PUBLIC API
makePartialLeastSquaresRegression <- function(xNamesVector, yNamesVector,
                                              XDataFrame, YDataFrame,
                                              treiningTestBoundary = 0.85, ncompValue = 10) {
    # Where XData = transcriptomicsData and YData = lipidomicsData.
    XData <- data.frame(XDataFrame[!duplicated(XDataFrame[1]), ], row.names = 1)
    YData <- data.frame(YDataFrame[!duplicated(YDataFrame[1]), ], row.names = 1)

    transposedXData <- as.data.frame(t(XData))
    transposedYData <- as.data.frame(t(YData))

    interX <- intersect(colnames(transposedXData), xNamesVector)
    interY <- intersect(colnames(transposedYData), yNamesVector)

    X <- transposedXData[as.character(interX)]
    Y <- transposedYData[as.character(interY)]

    Xmelt <- I(as.matrix(X))
    Ymelt <- I(as.matrix(Y))

    combined <- data.frame(X = I(Xmelt), Y = I(Ymelt))

    trainingRows <- ceiling(treiningTestBoundary * nrow(combined))
    combinedToTraining <- combined[1:trainingRows,]
    combinedToTest <- combined[(trainingRows + 1):nrow(combined),]

    PLSResults <- tryCatch(
        {
            PLSResultsFromMatrixInDF <- plsr(Y ~ X, data = combinedToTraining)

            PlsPredict <- predict(PLSResultsFromMatrixInDF, ncomp = ncompValue, newdata = combinedToTest)
            PlsTestRmsep <- RMSEP(PLSResultsFromMatrixInDF, newdata = combinedToTest)
            varianceExplained <- pls::explvar(PLSResultsFromMatrixInDF)

            #TODO: TRAINING and TEST is required!
            PLSResultsList <- list("training" = PLSResultsFromMatrixInDF,
                               "varianceExplained" = varianceExplained,
                               "test" = PlsPredict,
                               "testRmsep" = PlsTestRmsep)
            PLSResultsList
        },
        error = function(cond) {
            message("ONION - Included PLS can not solve task.")
            message("Original message (pls):")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning = function(cond) {
            message("ONION - Included PLS present warning.")
            message("Original message (pls):")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally = {
            message("ONION - PLS (pls) finished.")
        }
    )
    PLSResults
}

# NEW PUBLIC API
makePermutationTestOnPLS <- function(XDataFrame, YDataFrame, numberOfRowsForTestOnX, numberOfRowsForTestOnY, numberOfIterations = 100, countedPLS) {
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(numberOfRowsForTestOnX)
    print(numberOfRowsForTestOnY)
    # print()
    # print()
    vectorOfVarianceExplained <- as.numeric();
    for (i in 1:numberOfIterations) {
        xNV <- as.character(XDataFrame[sample(nrow(XDataFrame), numberOfRowsForTestOnX), ][,1])
        yNV <- as.character(YDataFrame[sample(nrow(YDataFrame), numberOfRowsForTestOnY), ][,1])
        plsResult <- ONION::makePartialLeastSquaresRegression(xNamesVector = xNV, yNamesVector = yNV, XDataFrame = XDataFrame, YDataFrame = YDataFrame)
        print("++++++++++++++++++++++++++++++++++")
        # print.default(ccaResult)
        if (is.na(plsResult) || is.null(plsResult)) {

        } else {
            vectorOfVarianceExplained <- c(vectorOfVarianceExplained, as.numeric(plsResult$varianceExplained[1]))
        }
    }

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    # print.default(countedCCA)
    print("*****************************")
    # print(countedCCA$xrd)
    # print(countedCCA$yrd)
    print(vectorOfVarianceExplained)
    print(mean(vectorOfVarianceExplained))
    print(countedPLS)

    meanOfVarianceExplained <- NA;
    if (0 != length(vectorOfVarianceExplained)) {
        meanOfVarianceExplained <- mean(vectorOfVarianceExplained);
    }

    countedPLSRecord <- NA;
    if (is.na(countedPLS) || is.null(countedPLS)) {

    } else {
        countedPLSRecord <- countedPLS$varianceExplained[1];
    }

    testResult <- list("countedPLSOnX" = countedPLSRecord,
                       "meanOnVarianceExplained" = meanOfVarianceExplained)
    testResult
}

# NEW PUBLIC API
makePLSOnGroups <- function(groupsDefinitionDF, mappingDF, leftMappingColumnName = 'root', rightMappingColumnName = 'genesSymbolsFromEnsemble', groupsDataDF, mappingDataDF){
    ddply(.data = groupsDefinitionDF['Molecules'], .(Molecules), .fun = function(dfElement) {
        print("???????????????????")
        print(dfElement)
        rightSideIdsToAnalys <- unlist(strsplit(as.character(dfElement), split = " "));
        print(rightSideIdsToAnalys)
        leftSideIdsToAnalys <- mappingDF[mappingDF[[leftMappingColumnName]] %in% rightSideIdsToAnalys,][[rightMappingColumnName]]
        leftSideIdsToAnalys <- unique(unlist(leftSideIdsToAnalys))

        plsResults <- makePartialLeastSquaresRegression(
            xNamesVector = leftSideIdsToAnalys,
            yNamesVector = rightSideIdsToAnalys,
            XDataFrame = mappingDataDF,
            YDataFrame = groupsDataDF)

        print("######################################")
        # print(ccaResults)
        print("---------------------------")
        # print.default(ccaResults)
        #TODO : Use user defined column name instead of symbol. ChEBI column too.
        numberOfRowsForTestOnX <- nrow(mappingDataDF[mappingDataDF$symbol %in% leftSideIdsToAnalys, ])
        numberOfRowsForTestOnY <- nrow(groupsDataDF[groupsDataDF$ChEBI %in% rightSideIdsToAnalys, ])
        parmutationTestResult <- makePermutationTestOnPLS(XDataFrame = mappingDataDF, YDataFrame = groupsDataDF,
                                                          numberOfRowsForTestOnX = numberOfRowsForTestOnX,
                                                          numberOfRowsForTestOnY = numberOfRowsForTestOnY,
                                                          numberOfIterations = 50, countedPLS = plsResults);

        dfWithCca <- data.frame('right' = I(list(rightSideIdsToAnalys)),
                                'left' = I(list(leftSideIdsToAnalys)),
                                'plsResults' = I(list(plsResults)),
                                'plsPermutationTestResults' = I(list(parmutationTestResult)))
        dfWithCca
    })
}

# NEW PUBLIC API
plotRmsepForPLS <- function(PLSResult) {
    plot(pls::RMSEP(PLSResult), legendpos = "topright")
}


# NEW PUBLIC API
plotRegression <- function(PLSResult, ncompValue) {
    plot(PLSResult, ncomp = ncompValue, asp = 1, line = TRUE)
}

# TODO: Analiza różnicowa, differencial analysis.

# TODO: Check exceptions hadling in methods.

# TODO: Check plots. :)
makePLSCharts <- function(PLS) {
    png("PLS_loadings.png", width = 640, height = 480)
        par(mfrow = c(2,2))
        biplot(PLS, which = "x") # Default
        biplot(PLS, which = "y")
        biplot(PLS, which = "scores")
        biplot(PLS, which = "loadings")
    dev.off()
}


# NEW PUBLIC API
createFunctionalInteractionsDataFrame <- function(chebiToReactomeDataFrame, singleIdColumnName = 'ontologyId', idsListColumnName = 'ensembleIds') {
    functionalInteractionsDataFrame <- ddply(.data = chebiIdsToReactomePathways, c(singleIdColumnName), .fun = function(dfElement) {
        functionalInteractionsRows <- adply(.data = dfElement[1,c(idsListColumnName)][[1]], .margins = 1, dfff = dfff, .fun = function(listElement, dfff) {
            functionalInteractionsRow <- data.frame("Gene1" = dfElement[1, c(singleIdColumnName)],
                               "Gene2" = listElement,
                               "Annotation" = "reactome",
                               "Direction" = "-",
                               "Score" = 1.00,
                               stringsAsFactors = FALSE)
            functionalInteractionsRow
        })
        functionalInteractionsRows
    })
    functionalInteractionsDataFrame[,c("Gene1", "Gene2", "Annotation", "Direction", "Score")]
}