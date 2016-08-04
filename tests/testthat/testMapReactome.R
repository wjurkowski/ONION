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
    #Filtering, analiza składowych, zmiana kolumn DF w baseData.
    clusteredSmallMolecules <- ONION::clusterUsingOntology(chebiIdsDataFrame = baseData,
                                ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology)
    head(clusteredSmallMolecules)

    mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(clusteredSmallMolecules)
    head(mergedSmallMolecules)

    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules, organismTaxonomyId = '9606')
    head(chebiIdsToReactomePathways)

    functionalInteractions <- ONION::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)
    head(functionalInteractions)

    chebiIdsToReactomePathwaysSmall <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[c("ontologyId")], organismTaxonomyId = '9606')
    head(chebiIdsToReactomePathwaysSmall)

    chebiIdsToReactomePathwaysAndToStringNeighbours <- ONION::getStringNeighbours(chebiIdsToReactomePathways)
    head(chebiIdsToReactomePathwaysAndToStringNeighbours)

    gmtGroupsFilePath <- paste(find.package("ONION"),"/example/nm-groups.txt", sep = "")
    groups <- ONION::groupUsingUserDefinition(pathToFileWithGroupDefinition = gmtGroupsFilePath)

    # JOIN on input data and Reactome result, String results.
    # INFO : Be careful Ewa with 73705 -> 16015. Look at clusteredSmallMolecules and mergedSmallMolecules.
    lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == 27432,]$root
    lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == 73705,]$root
    joinLip <- c(lip1, lip2)

    reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$chebiId == 27432,]$gensSymbols[[1]]
    reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$chebiId == 16015,]$gensSymbols[[1]]
    joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

    stringTrans1 <- chebiIdsToReactomePathwaysAndToStringNeighbours[
            chebiIdsToReactomePathwaysAndToStringNeighbours$chebiId == 27432,
        ]$stringGensSymbols[[1]]
    stringTrans2 <- chebiIdsToReactomePathwaysAndToStringNeighbours[
            chebiIdsToReactomePathwaysAndToStringNeighbours$chebiId == 16015,
        ]$stringGensSymbols[[1]]
    joinStringTrans <- c(stringTrans1, stringTrans2)[!duplicated(c(stringTrans1, stringTrans2))]


    pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
    pathToExampleFileWithYData <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")

    XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
    YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

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

    plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults1)

    plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults2)


    PLSResult1 <- ONION::makePartialLeastSquaresRegression(
        joinRecatomeTrans,
        joinLip,
        XDataFrame = XDF,
        YDataFrame = YDF)

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

