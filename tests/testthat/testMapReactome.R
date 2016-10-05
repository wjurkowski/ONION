#Set working directory. Solve relative path problem in tests.
#TODO Add different behavior on CHECK and TEST.
#setwd(paste(getwd(), "/.."))
#print(getwd())

#source("../R/mapReactome.R", chdir = TRUE)


set_up <- function() {
    source("R/mapReactome.R")
}

tear_down <- function() {
    cat("\014")
}

test_that("NewOnionApiWorkflow test", {
    #given
    print('*********given*********')
    # Basic flow. Chebi -> Reactome -> String.
    pathToFileWithChebiIds <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")
    baseData <- read.table(pathToFileWithChebiIds, header = TRUE)
    head(baseData)
    #Filtering, analiza składowych, zmiana kolumn DF w baseData.
    clusteredSmallMolecules <- ONION::clusterUsingOntology(chebiIdsDataFrame = baseData,
                                                           rootColumnName = "ChEBI",
                                                           ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology)
    head(clusteredSmallMolecules)

    mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(clusteredSmallMolecules, rootColumnName = 'root')
    head(mergedSmallMolecules)

    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules, organismTaxonomyId = '9606')
    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[3,], organismTaxonomyId = '9606')
    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[4,], organismTaxonomyId = '9606')
    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[1,], organismTaxonomyId = '9606')
    head(chebiIdsToReactomePathways)

    functionalInteractions <- ONION::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)
    head(functionalInteractions)

    chebiIdsToReactomePathwaysSmall <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[c("ontologyId")], organismTaxonomyId = '9606', rootColumnName = NULL)
    head(chebiIdsToReactomePathwaysSmall)


    chebiIdsToReactomePathwaysAndToStringNeighbours <- ONION::getStringNeighbours(chebiIdsToReactomePathways[1,],
                                                                                  stringOrganismId = 9606,
                                                                                  stringDbVersion = "10",
                                                                                  idsColumnName = 'ontologyId',
                                                                                  rootColumnName = NULL,
                                                                                  listOfEnsembleIdColumnName = 'ensembleIds')
    head(chebiIdsToReactomePathwaysAndToStringNeighbours)

    gmtGroupsFilePath <- paste(find.package("ONION"),"/example/nm-groups.txt", sep = "")
    groups <- ONION::groupUsingUserDefinition(pathToFileWithGroupDefinition = gmtGroupsFilePath)

    # JOIN on input data and Reactome result, String results.
    # INFO : Be careful Ewa with 73705 -> 16015. Look at clusteredSmallMolecules and mergedSmallMolecules.
    lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:27432",]$root
    lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
    joinLip <- c(as.character(lip1), as.character(lip2))

    reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:27432",]$genesSymbols[[1]]
    reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:16015",]$genesSymbols[[1]]
    joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

    stringTrans1 <- chebiIdsToReactomePathwaysAndToStringNeighbours[
            chebiIdsToReactomePathwaysAndToStringNeighbours$chebiId == 27432,
        ]$stringGenesSymbols[[1]]
    stringTrans2 <- chebiIdsToReactomePathwaysAndToStringNeighbours[
            chebiIdsToReactomePathwaysAndToStringNeighbours$chebiId == 16015,
        ]$stringGenesSymbols[[1]]
    joinStringTrans <- c(stringTrans1, stringTrans2)[!duplicated(c(stringTrans1, stringTrans2))]


    pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
    pathToExampleFileWithYData <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")

    XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
    YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

    dim(YDF)
    dim(XDF)

    ccaResults1 <- ONION::makeCanonicalCorrelationAnalysis(
        xNamesVector = joinRecatomeTrans,
        yNamesVector = joinLip,
            XDataFrame = XDF,
            YDataFrame = YDF)

    ccaResults2 <- ONION::makeCanonicalCorrelationAnalysis(
        joinStringTrans,
        joinLip,
        XDataFrame = XDF,
        YDataFrame = YDF)

    groups['Molecules']

    rooot <- 'root'
    chebiIdsToReactomePathways[rooot]
    chebiIdsToReactomePathways[[rooot]]
    chebiIdsToReactomePathways$root

    mccReactome <- ONION::makeCCAOnGroups(groupsDefinitionDF = groups, mappingDF = chebiIdsToReactomePathways, groupsDataDF = YDF, mappingDataDF = XDF)

    permutationTestsResults <- ONION::makePermutationTestOnCCA(XDataFrame = XDF, YDataFrame = YDF, 17, 2, 50, countedCCA = ccaResults1)

    plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome[4,]$ccaResults[[1]])

    mccReactome[4,]$ccaPermutationTestResults

    p <- recordPlot()

    plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome[4,]$ccaResults[[1]])
    p

    plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome[1,]$ccaResults[[1]])

    mccString <- ONION::makeCCAOnGroups(groupsDefinitionDF = groups, mappingDF = chebiIdsToReactomePathwaysAndToStringNeighbours, groupsDataDF = YDF, mappingDataDF = XDF, rightMappingColumnName = 'stringGenesSymbols')
    mccString[7,]$ccaResults

    plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults1)

    plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults2)


    PLSResult1 <- NULL

    PLSResult1 <- ONION::makePartialLeastSquaresRegression(
        joinRecatomeTrans,
        joinLip,
        XDataFrame = XDF,
        YDataFrame = YDF)

    groupPlsReactome <- ONION::makePLSOnGroups(groupsDefinitionDF = groups, mappingDF = chebiIdsToReactomePathways, groupsDataDF = YDF, mappingDataDF = XDF)

    as.numeric(PLSResult1$varianceExplained[1])

    explanatoryVariance <- pls::explvar(PLSResult1$training)


    PLSResult2 <- ONION::makePartialLeastSquaresRegression(
        joinStringTrans,
        joinLip,
        XDataFrame = XDF,
        YDataFrame = YDF)



    summary(PLSResult1$training)
    summary(PLSResult2$training)

    ONION::plotRmsepForPLS(PLSResult1$training)
    ONION::plotRmsepForPLS(PLSResult2$training)

    ONION::plotRegression(PLSResult1$training, ncompValue = 10)
    ONION::plotRegression(PLSResult2$training, ncompValue = 10)

    PLSResult1$test
    PLSResult1$testRmsep

    PLSResult2$test
    PLSResult2$testRmsep

    # IMPORTANT
    # Counterexample of STRINGdb.
    stringOrganismId = 9606
    stringDbVersion = "10"
    library(STRINGdb)
    string_db <- STRINGdb$new( version = stringDbVersion, species = stringOrganismId)
    stringId1 <- string_db$mp(c("ENSG00000115263"))
    returnNeighbourVector <- string_db$get_neighbors(stringId1)
    duplicated(returnNeighbourVector)
})

