
OmicsON::setUpReactomeMapping(ChEBI2ReactomeFileURL = "https://reactome.org/download/current/ChEBI2Reactome.txt",
                              Ensembl2ReactomeFileURL = "https://reactome.org/download/current/Ensembl2Reactome.txt",
                              UniProt2ReactomeFileURL = "https://reactome.org/download/current/UniProt2Reactome.txt")

pathToFileWithLipidomicsData <- system.file(package="OmicsON", "extdata", "nm-lipidomics.txt")

lipidomicsInputData <- read.table(pathToFileWithLipidomicsData, header = TRUE)
lipidomicsInputDf <- head(lipidomicsInputData, 6)
knitr::kable(lipidomicsInputDf[1:7], caption = "Lipidomisc data")

pathToFileWithTranscriptomicsData <- system.file(package="OmicsON", "extdata", "nm-transcriptomics.txt")

transcriptomicsInputData <- read.table(pathToFileWithTranscriptomicsData, header = TRUE)
transcriptomicsInputDf <- head(transcriptomicsInputData, 6)
knitr::kable(transcriptomicsInputDf[1:7], caption = "Transcriptomics data")

decReac <- OmicsON::decorateByReactomeData(chebiMoleculesDf = lipidomicsInputData,
                                chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')

knitr::kable(lipidomicsInputData[c(2, 3, 6),])
knitr::kable(decReac[c(2, 12, 18),])
for(x in decReac[3]){
    ddply(decReac, .(ontologyId), function(x) {
        print(length(unlist(x$ensembleIds)))
    })
    print(length(x))
}
length(decReac[3])


# Easy API
# DONE : root column jako pierwsza.
# DONE : ontoloogyId - dac te same id co w root ale puste mapowania.
# FIXED : BUG : Chebi id function - sprawdzić czy działa.


hist(unlist(lipidomicsInputData[1:20,-c(1)]*100), breaks=350)
hist(AirPassengers)

dataDecoratedByReactome <- decorateByReactomeData(chebiMoleculesDf = lipidomicsInputData,
                                                  chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')
dataDecoratedByReactome

chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]

decorateByStringData <- function(chebiMoleculesDecoratedByReactomeDf, organismTaxonomyId = '9606') {
    chebiIdsToReactomePathwaysAndToStringNeighbours <- OmicsON::getStringNeighbours(
        chebiIdsToReactomePathways = chebiMoleculesDecoratedByReactomeDf,
        stringOrganismId = as.numeric(organismTaxonomyId),
        stringDbVersion = "10",
        idsColumnName = 'ontologyId',
        rootColumnName = NULL,
        listOfEnsembleIdColumnName = 'ensembleIds')
}

library(STRINGdb)
string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=0, input_directory="" )

source("https://bioconductor.org/biocLite.R")
biocLite("STRINGdb")

dataDecoratedByString <- decorateByStringData(
    chebiMoleculesDecoratedByReactomeDf = decReac[decReac$ontologyId == "CHEBI:15756",]
)

dataDecoratedByString["stringIds"]$stringIds
dataDecoratedByString["ensembleIds"]$ensembleIds


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


#select small molecules
decReac[c(12),]
lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:15756",]$root
lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
joinLip <- c(as.character(lip1), as.character(lip2))

lipidomics <- c("CHEBI:73705", "CHEBI:15756")
#use Reactome genes mapped to selected small molecules
reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:28875",]$genesSymbolsFromEnsemble[[1]]
joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

#functionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)
functionalInteractions <- OmicsON::createFunctionalInteractionsDataFrame(dataDecoratedByReactome)

pathToExampleFileWithXData <- paste(find.package("OmicsON"),"/example/nm-transcriptomics.txt", sep = "")
pathToExampleFileWithYData <- paste(find.package("OmicsON"),"/example/nm-lipidomics.txt", sep = "")


#NEW API
decReac

transcriptomics <- decReac[decReac$root == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
transcriptomics <- decReac[decReac$root == "CHEBI:72850",]$genesSymbolsFromEnsemble[[1]]
transcriptomics <- decReac[decReac$root == "CHEBI:28875",]$genesSymbolsFromEnsemble[[1]]
lipidomics <- c("CHEBI:72850", "CHEBI:15756", "CHEBI:28875")

XDataFrame <- transcriptomicsInputData
YDataFrame <- lipidomicsInputData

XData <- data.frame(XDataFrame[!duplicated(XDataFrame[1]), ], row.names = 1)
YData <- data.frame(YDataFrame[!duplicated(YDataFrame[1]), ], row.names = 1)

transposedXData <- as.data.frame(t(XData))
transposedYData <- as.data.frame(t(YData))

xNamesVector <- transcriptomics
yNamesVector <- lipidomics
interX <- intersect(colnames(transposedXData), xNamesVector)
interY <- intersect(colnames(transposedYData), yNamesVector)

X <- transposedXData[as.character(interX)]
Y <- transposedYData[as.character(interY)]

#BUG : nrow() not length!!!
length(X)
length(Y)
ccaResults1 <- yacca::cca(X, Y)
ccaResults1 <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = transcriptomics,
    yNamesVector = lipidomics,
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)

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
str(PLSResult1)

groupPlsReactome <- OmicsON::makePLSOnGroups(
    groupsDefinitionDF = groups,
    mappingDF = chebiIdsToReactomePathwaysWithRoot,
    groupsDataDF = YDF,
    mappingDataDF = XDF)

OmicsON::plotRmsepForPLS(groupPlsReactome$plsResults[[1]]$training)

OmicsON::plotRegression(groupPlsReactome$plsResults[[1]]$training, ncompValue = 10)


