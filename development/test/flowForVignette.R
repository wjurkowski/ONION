
pathToFileWithLipidomicsData <- system.file(package="OmicsON", "example", "nm-lipidomics.txt")

lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
lipidomicsInputDf <- head(lipidomicsInputData, 6)
knitr::kable(lipidomicsInputDf[1:7], caption = "Lipidomisc data")

pathToFileWithTranscriptomicsData <- system.file(package="OmicsON", "example", "nm-transcriptomics.txt")

transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)
transcriptomicsInputDf <- head(transcriptomicsInputData, 6)
knitr::kable(transcriptomicsInputDf[1:7], caption = "Transcriptomics data")

clusteredSmallMolecules <- OmicsON::clusterUsingOntology(
    chebiIdsDataFrame = lipidomicsInputDf,
    rootColumnName = "ChEBI",
    ontologyRepresentatnion = OmicsON::firstExistsInReactomeChebiOntology)
knitr::kable(head(clusteredSmallMolecules, 6))

mergedSmallMolecules <- OmicsON::mergeChEBIOntologyWithChildFavoring(
    clusteredSmallMolecules,
    rootColumnName = 'root')

knitr::kable(head(mergedSmallMolecules, 6))

knitr::kable(data.frame(
    mergedSmallMolecules[1:4, c("ontologyId"), drop = FALSE],
    XYZ = c("1","2","3","4"))
)

chebiIdsToReactomePathways <- OmicsON::mapReactomePathwaysUnderOrganism(
    chebiOntologyIds = mergedSmallMolecules[, c("ontologyId"), drop = FALSE],
    organismTaxonomyId = '9606',
    idsColumnName = "ontologyId",
    rootColumnName = NULL)
chebiIdsToReactomePathwaysWithRoot <- OmicsON::mapReactomePathwaysUnderOrganism(
    chebiOntologyIds = mergedSmallMolecules[, c("ontologyId", "root"), drop = FALSE],
    organismTaxonomyId = '9606',
    idsColumnName = "ontologyId",
    rootColumnName = "root")


oneRowDf <- chebiIdsToReactomePathways[2,]
rownames(oneRowDf) <- NULL
knitr::kable(oneRowDf)


chebiIdsToReactomePathwaysAndToStringNeighbours <- OmicsON::getStringNeighbours(
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


gmtGroupsFilePath <- system.file(package="OmicsON", "example", "nm-groups.txt")

groups <- OmicsON::readGroupsAsDf(pathToFileWithGroupDefinition = gmtGroupsFilePath)
knitr::kable(groups)

#select small molecules
lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:15756",]$root
lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
joinLip <- c(as.character(lip1), as.character(lip2))

#use Reactome genes mapped to selected small molecules
reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:28875",]$genesSymbolsFromEnsemble[[1]]
joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

functionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)


pathToExampleFileWithXData <- paste(find.package("OmicsON"),"/example/nm-transcriptomics.txt", sep = "")
pathToExampleFileWithYData <- paste(find.package("OmicsON"),"/example/nm-lipidomics.txt", sep = "")

XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

ccaResults1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = joinRecatomeTrans,
    yNamesVector = joinLip,
    XDataFrame = XDF,
    YDataFrame = YDF)

OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResults1)

mccReactome <- OmicsON::makeCCAOnGroups(
    groupsDefinitionDF = groups,
    mappingDF = chebiIdsToReactomePathwaysWithRoot,
    groupsDataDF = YDF,
    mappingDataDF = XDF)

OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = mccReactome$ccaResults[[1]])

PLSResult1 <- OmicsON::makePartialLeastSquaresRegression(
    joinRecatomeTrans,
    joinLip,
    XDataFrame = XDF,
    YDataFrame = YDF)

OmicsON::plotRmsepForPLS(PLSResult1$training)

OmicsON::plotRegression(PLSResult1$training, ncompValue = 10)

groupPlsReactome <- OmicsON::makePLSOnGroups(
    groupsDefinitionDF = groups,
    mappingDF = chebiIdsToReactomePathwaysWithRoot,
    groupsDataDF = YDF,
    mappingDataDF = XDF)

OmicsON::plotRmsepForPLS(groupPlsReactome$plsResults[[1]]$training)

OmicsON::plotRegression(groupPlsReactome$plsResults[[1]]$training, ncompValue = 10)
