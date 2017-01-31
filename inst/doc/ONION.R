## ----global_options, include=TRUE----------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, error = TRUE)

## ---- results = 'asis'---------------------------------------------------
    pathToFileWithLipidomicsData <- paste(
        find.package("ONION"),
        "/example/nm-lipidomics.txt", 
        sep = "")
    lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
    lipidomicsInputDf <- head(lipidomicsInputData, 6)
    knitr::kable(lipidomicsInputDf[1:7], caption = "Lipidomisc data")
    
    pathToFileWithTranscriptomicsData <- paste(
        find.package("ONION"),
        "/example/nm-transcriptomics.txt", 
        sep = "")
    transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)
    transcriptomicsInputDf <- head(transcriptomicsInputData, 6)
    knitr::kable(transcriptomicsInputDf[1:7], caption = "Transcriptomics data")

## ---- echo=TRUE, results='asis'------------------------------------------
    clusteredSmallMolecules <- ONION::clusterUsingOntology(
        chebiIdsDataFrame = lipidomicsInputDf,
        rootColumnName = "ChEBI",
        ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology)

## ---- echo=FALSE, results='asis'-----------------------------------------
    knitr::kable(head(clusteredSmallMolecules, 6))

## ---- echo=TRUE, results='asis'------------------------------------------
    mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(
        clusteredSmallMolecules, 
        rootColumnName = 'root')

## ---- echo=FALSE, results='asis'-----------------------------------------
    knitr::kable(head(mergedSmallMolecules, 6))

## ---- echo=FALSE, results='asis'-----------------------------------------
    knitr::kable(data.frame(
        mergedSmallMolecules[1:4, c("ontologyId"), drop = FALSE], 
        XYZ = c("1","2","3","4"))
    )

## ---- echo=TRUE, results='hide'------------------------------------------
    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(
        chebiOntologyIds = mergedSmallMolecules[, c("ontologyId"), drop = FALSE], 
        organismTaxonomyId = '9606', 
        idsColumnName = "ontologyId", 
        rootColumnName = NULL)
    chebiIdsToReactomePathwaysWithRoot <- ONION::mapReactomePathwaysUnderOrganism(
        chebiOntologyIds = mergedSmallMolecules[, c("ontologyId", "root"), drop = FALSE], 
        organismTaxonomyId = '9606', 
        idsColumnName = "ontologyId", 
        rootColumnName = "root")

## ---- echo=FALSE, results='asis'-----------------------------------------
    oneRowDf <- chebiIdsToReactomePathways[6,]
    rownames(oneRowDf) <- NULL
    knitr::kable(oneRowDf)

## ---- echo=TRUE, results='hide'------------------------------------------
chebiIdsToReactomePathwaysAndToStringNeighbours <- ONION::getStringNeighbours(
    chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",],
    stringOrganismId = 9606,
    stringDbVersion = "10",
    idsColumnName = 'ontologyId',
    rootColumnName = NULL,
    listOfEnsembleIdColumnName = 'ensembleIds')

## ---- echo=FALSE, results='asis'-----------------------------------------
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringIds[[1]] <- chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringIds[[1]][1:50]
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringGenesSymbols[[1]] <-
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringGenesSymbols[[1]][1:45]
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$ensembleIds[[1]] <-
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$ensembleIds[[1]][1:11]
    knitr::kable(
        chebiIdsToReactomePathwaysAndToStringNeighbours[
            chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]
    )

## ---- echo=TRUE, results='asis'------------------------------------------
    gmtGroupsFilePath <- paste(find.package("ONION"),"/example/nm-groups.txt", sep = "")
    groups <- ONION::readGroupsAsDf(pathToFileWithGroupDefinition = gmtGroupsFilePath)

## ---- echo=FALSE, results='asis'-----------------------------------------
    knitr::kable(groups)

## ---- echo=TRUE, results='asis'------------------------------------------
#select small molecules
lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:27432",]$root
lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
joinLip <- c(as.character(lip1), as.character(lip2))

#use Reactome genes mapped to selected small molecules
reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:16015",]$genesSymbolsFromEnsemble[[1]]
joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

## ---- echo=TRUE, results='asis'------------------------------------------
    functionalInteractions <- ONION::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)

## ---- echo=FALSE, results='asis'-----------------------------------------
    knitr::kable(head(functionalInteractions, 6))

## ---- echo=TRUE, results='asis'------------------------------------------
    pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
    pathToExampleFileWithYData <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")

    XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
    YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

    ccaResults1 <- ONION::makeCanonicalCorrelationAnalysis(
        xNamesVector = joinRecatomeTrans,
        yNamesVector = joinLip,
            XDataFrame = XDF,
            YDataFrame = YDF)

## ---- fig.show='hold', fig.width=6, fig.height=6-------------------------
    ONION::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults1)

## ---- echo=TRUE, results='hide'------------------------------------------
    mccReactome <- ONION::makeCCAOnGroups(
        groupsDefinitionDF = groups, 
        mappingDF = chebiIdsToReactomePathwaysWithRoot, 
        groupsDataDF = YDF, 
        mappingDataDF = XDF)


## ---- echo=FALSE, results='hide'-----------------------------------------
    permutationTestsResults <- ONION::makePermutationTestOnCCA(
        XDataFrame = XDF, 
        YDataFrame = YDF, 17, 2, 50, 
        countedCCA = ccaResults1)

## ---- echo=FALSE, results='hide'-----------------------------------------
mccReactome$Molecules[1]
mccReactome$right[[1]]
mccReactome$left[[1]]
mccReactome$ccaResults[[1]]
mccReactome$ccaPermutationTestResults[[1]]

## ---- fig.show='hold', fig.width=6, fig.height=6-------------------------
    ONION::plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome$ccaResults[[1]])

## ---- echo=FALSE, results='hide'-----------------------------------------
    PLSResult1 <- ONION::makePartialLeastSquaresRegression(
        joinRecatomeTrans,
        joinLip,
        XDataFrame = XDF,
        YDataFrame = YDF)

## ---- fig.show='hold', fig.width=6, fig.height=6-------------------------
    ONION::plotRmsepForPLS(PLSResult1$training)

## ---- fig.show='hold', fig.width=6, fig.height=6-------------------------
    ONION::plotRegression(PLSResult1$training, ncompValue = 10)

## ---- echo=TRUE, results='asis'------------------------------------------
    groupPlsReactome <- ONION::makePLSOnGroups(
        groupsDefinitionDF = groups, 
        mappingDF = chebiIdsToReactomePathwaysWithRoot, 
        groupsDataDF = YDF, 
        mappingDataDF = XDF)


## ---- echo=FALSE, results='hide'-----------------------------------------
groupPlsReactome$Molecules[1]
groupPlsReactome$right[[1]]
groupPlsReactome$left[[1]]
groupPlsReactome$plsResults[[1]]
groupPlsReactome$plsPermutationTestResults[[1]]

## ---- fig.show='hold', fig.width=6, fig.height=6, echo=TRUE--------------
    ONION::plotRmsepForPLS(groupPlsReactome$plsResults[[1]]$training)

## ---- fig.show='hold', fig.width=6, fig.height=6, echo=TRUE--------------
    ONION::plotRegression(groupPlsReactome$plsResults[[1]]$training, ncompValue = 10)

