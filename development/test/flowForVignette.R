
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

lipidomicsInputData[c(2, 3, 6),]
decReac[c(2, 12, 18),]


# Easy API
# DONE : root column jako pierwsza.
# DONE : ontoloogyId - dac te same id co w root ale puste mapowania.
# FIXED : BUG : Chebi id function - sprawdzić czy działa.

decSrtDb <- OmicsON::decorateByStringDbData(chebiIdsToReactomePathways = decReac, listOfEnsembleIdColumnName = 'ensembleIds')

decSrtDbUniProt <- OmicsON::decorateByStringDbData(chebiIdsToReactomePathways = decReac, listOfEnsembleIdColumnName = 'uniProtIds')


decSrtDbUniProtNameTest <- OmicsON::decorateByStringDbData(chebiIdsToReactomePathways = decReac[1,], listOfEnsembleIdColumnName = 'uniProtIds')


pathToExampleFileWithXData <- paste(find.package("OmicsON"),"/example/nm-transcriptomics.txt", sep = "")
pathToExampleFileWithYData <- paste(find.package("OmicsON"),"/example/nm-lipidomics.txt", sep = "")

XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
YDF <- read.table(pathToExampleFileWithYData, header = TRUE);

transcriptomicsInputData
lipidomicsInputData

typedTransData <- decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:28875"),"stringGenesSymbolsExpand"]
typedTransData <- decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:73705"),"stringGenesSymbolsExpand"]

typedTransData <- decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:28875"),"stringGenesSymbolsNarrow"]
typedTransData <- decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:73705"),"stringGenesSymbolsNarrow"]


trandData <- transcriptomicsInputData$symbol
intersect(typedTransData[[1]], as.character(trandData))
str(typedTransData)
str(trandData)

plyr::ddply(.data = decSrtDb, .variables = c("root"), .fun = function(dfRow) {
    data.frame("common" = I(list(intersect(dfRow[,"stringGenesSymbolsExpand"][[1]], as.character(trandData)))))
})

ccaResultsExpand <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:73705"),"stringGenesSymbolsExpand"][[1]],
    yNamesVector = c("CHEBI:73705"),
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)

ccaResultsNarrow <- OmicsON::makeCanonicalCorrelationAnalysis(
    xNamesVector = decSrtDb[decSrtDb[,"root"] %in% c("CHEBI:73705"),"stringGenesSymbolsNarrow"][[1]],
    yNamesVector = c("CHEBI:73705"),
    XDataFrame = transcriptomicsInputData,
    YDataFrame = lipidomicsInputData)



OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsExpand)
OmicsON::plotCanonicalCorrelationAnalysisResults(ccaResults = ccaResultsNarrow)








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


