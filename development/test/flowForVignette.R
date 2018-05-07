
pathToFileWithLipidomicsData <- system.file(package="ONION", "example", "nm-lipidomics.txt")

lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
lipidomicsInputDf <- head(lipidomicsInputData, 6)
knitr::kable(lipidomicsInputDf[1:7], caption = "Lipidomisc data")

pathToFileWithTranscriptomicsData <- system.file(package="ONION", "example", "nm-transcriptomics.txt")

transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)
transcriptomicsInputDf <- head(transcriptomicsInputData, 6)
knitr::kable(transcriptomicsInputDf[1:7], caption = "Transcriptomics data")

clusteredSmallMolecules <- ONION::clusterUsingOntology(
    chebiIdsDataFrame = lipidomicsInputDf,
    rootColumnName = "ChEBI",
    ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology)
knitr::kable(head(clusteredSmallMolecules, 6))

mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(
    clusteredSmallMolecules,
    rootColumnName = 'root')

knitr::kable(head(mergedSmallMolecules, 6))

knitr::kable(data.frame(
    mergedSmallMolecules[1:4, c("ontologyId"), drop = FALSE],
    XYZ = c("1","2","3","4"))
)

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


oneRowDf <- chebiIdsToReactomePathways[2,]
rownames(oneRowDf) <- NULL
knitr::kable(oneRowDf)


chebiIdsToReactomePathwaysAndToStringNeighbours <- ONION::getStringNeighbours(
    chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",],
    stringOrganismId = 9606,
    stringDbVersion = "10",
    idsColumnName = 'ontologyId',
    rootColumnName = NULL,
    listOfEnsembleIdColumnName = 'ensembleIds')

chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringIds[[1]] <-
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringIds[[1]][1:50]
chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringGenesSymbols[[1]] <-
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$stringGenesSymbols[[1]][1:45]
chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$ensembleIds[[1]] <-
    chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]$ensembleIds[[1]][1:11]
knitr::kable(
    chebiIdsToReactomePathwaysAndToStringNeighbours[
        chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",]
)


gmtGroupsFilePath <- system.file(package="ONION", "example", "nm-groups.txt")
# gmtGroupsFilePath <- paste(find.package("ONION"),"/example/nm-groups.txt", sep = "")
groups <- ONION::readGroupsAsDf(pathToFileWithGroupDefinition = gmtGroupsFilePath)
knitr::kable(groups)

#select small molecules
lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:15756",]$root
lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
joinLip <- c(as.character(lip1), as.character(lip2))

#use Reactome genes mapped to selected small molecules
reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:28875",]$genesSymbolsFromEnsemble[[1]]
joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

functionalInteractions <- ONION::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)


pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
pathToExampleFileWithYData <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")

XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

ccaResults1 <- ONION::makeCanonicalCorrelationAnalysis(
    xNamesVector = joinRecatomeTrans,
    yNamesVector = joinLip,
    XDataFrame = XDF,
    YDataFrame = YDF)

ONION::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults1)

mccReactome <- ONION::makeCCAOnGroups(
    groupsDefinitionDF = groups,
    mappingDF = chebiIdsToReactomePathwaysWithRoot,
    groupsDataDF = YDF,
    mappingDataDF = XDF)

ONION::plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome$ccaResults[[1]])

PLSResult1 <- ONION::makePartialLeastSquaresRegression(
    joinRecatomeTrans,
    joinLip,
    XDataFrame = XDF,
    YDataFrame = YDF)

ONION::plotRmsepForPLS(PLSResult1$training)

ONION::plotRegression(PLSResult1$training, ncompValue = 10)

groupPlsReactome <- ONION::makePLSOnGroups(
    groupsDefinitionDF = groups,
    mappingDF = chebiIdsToReactomePathwaysWithRoot,
    groupsDataDF = YDF,
    mappingDataDF = XDF)

ONION::plotRmsepForPLS(groupPlsReactome$plsResults[[1]]$training)

ONION::plotRegression(groupPlsReactome$plsResults[[1]]$training, ncompValue = 10)
