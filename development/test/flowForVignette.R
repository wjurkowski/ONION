
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

library(ONION)

decReac <- OmicsON::decorateByReactomeData(chebiMoleculesDf = lipidomicsInputData,
                                chebiIdsColumnName = "ChEBI", organismTaxonomyId = '9606')

knitr::kable(lipidomicsInputData[c(2, 3, 6),])
knitr::kable(decReac[c(2, 12, 18),])



# Easy API
# DONE : root column jako pierwsza.
# DONE : ontoloogyId - dac te same id co w root ale puste mapowania.
# FIXED : BUG : Chebi id function - sprawdzić czy działa.
library(STRINGdb)

chebiIdsToReactomePathways <- decReac
stringOrganismId = as.numeric('9606')
stringDbVersion = "10"
idsColumnName = 'ontologyId'
rootColumnName = 'root'
listOfEnsembleIdColumnName = 'ensembleIds'
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(rootColumnName, idsColumnName)
    }

    string_db <- STRINGdb$new(
        version = stringDbVersion,
        species = stringOrganismId,
        input_directory = path.expand("~"))
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[!chebiIdsToReactomePathways[idsColumnName] == '', ]
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[c(1,19),]
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways
    dfElement <- chebiIdsToRealReactomePathways[1,]
    dfElement <- chebiIdsToRealReactomePathways[19,]
    dfWithString <- ddply(.data = chebiIdsToRealReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
        extendedByStringAsVector <- character(length = 0)
        stringGenesSymbols <- character(length = 0)
        if (0 == length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
        } else {
            # diff_exp_example2 <- data.frame("translate" = decReac[,"ensembleIds"][[1]])
            # example2_mapped <- string_db$map( diff_exp_example2, "translate", removeUnmappedRows = TRUE )
            # example2_mapped[,"STRING_id"]

            toTranslate <- data.frame("translate" = dfElement[1,listOfEnsembleIdColumnName][[1]])
            translated <- string_db$map( toTranslate, "translate", removeUnmappedRows = TRUE )
            stringId1 <- translated[,"STRING_id"]
# cat("\014")
            # proteinIds <- getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(
            #     dfElement[1,listOfEnsembleIdColumnName], organismTaxonomyId = stringOrganismId
            # )
            # stringId1 <- string_db$mp(proteinIds)
            #New Approach!!!
            stringGraph <- string_db$get_graph()


            extendedByString <- plyr::ddply(.data = translated, .variables = c("STRING_id"), .fun = function(r) {
                data.frame("res" = I(list(igraph::neighbors(stringGraph, r[,"STRING_id"])$name)))
            })

            extendedByStringAsVector <- unique(unlist(extendedByString[,"res"]))

            # interSect <- unique(ssssaaaa[,"res"][[1]])
            # for(i in length(ssssaaaa[,"res"])) {
            #     interSect <- intersect(interSect, ssssaaaa[,"res"][[i]])
            # }


            # returnNeighbourVector <- string_db$get_neighbors(stringId1)
            ensembleIdsFromStringDb <- mapFromStringIdsToEnsembleIds(extendedByStringAsVector)
            stringGenesSymbols <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDb, organismTaxonomyId = stringOrganismId)
        }
        dffff <- data.frame('ensembleIds' = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(extendedByStringAsVector))),
                            'stringGenesSymbols' = I(list(unique(stringGenesSymbols))) )
        dffff
    })
    dfWithString
}






















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
data(diff_exp_example1)
head(diff_exp_example1)
diff_exp_example2 <- data.frame("gene" = dfElement[,"genesSymbolsFromEnsemble"][[1]])

diff_exp_example2 <- data.frame("gene" = dfElement[,"ensembleIds"][[1]])

example2_mapped <- string_db$map( diff_exp_example2, "gene", removeUnmappedRows = TRUE )
hits <- example2_mapped$STRING_id[1:400]
getOption("SweaveHooks")[["fig"]]()
string_db$plot_network( hits )

tp53 = string_db$mp( "tp53" )
atm = string_db$mp( "atm" )
neigh <- string_db$get_neighbors( c(tp53, atm) )
neigh

string_db <- STRINGdb$new( version="10", species=9606,
                           score_threshold=0, input_directory="" )
hits <- example2_mapped$STRING_id[3]
neigh <- "dupa"
neigh <- string_db$get_neighbors(string_ids = "9606.ENSP00000003084")
neigh
interr <- string_db$get_interactions(string_ids = "9606.ENSP00000005178")
interr$from
interr$to

string_proteins <- string_db$get_proteins()
ggg <- string_db$get_graph()
string_proteins
str(ggg)

cat("\014")
ggg[1,]
ggg[2,]

vvvv<-V(ggg)
asas1 <- igraph::neighbors(ggg, "9606.ENSP00000003084")
asas2 <- igraph::neighbors(ggg, "9606.ENSP00000005178")
unlist(asas1$name)
asas2$name
intersect(asas1$name, asas2$name)
asas3 <- igraph::neighbors(ggg, c("9606.ENSP00000005178","9606.ENSP00000003084"))
asas3 <- igraph::neighbors(ggg, c("9606.ENSP00000003084","9606.ENSP00000005178"))
asas1
asas2
asas3$name
unlist(asas2)
str(asas2)
asas1$name
E(ggg)



library(igraph)



colnames(interr) <- c("proteinA", "proteinB")
interScore <- cbind(interr, "score" = c(1.0))[c(1,2,17)]
interactions_benchmark = string_db$benchmark_ppi(
    interScore, pathwayType = "KEGG",
    max_homology_bitscore = 60, precision_window = 400, exclude_pathways = "blacklist")


data(interactions_example)
head(interactions_example)
interactions_benchmark = string_db$benchmark_ppi(
    interactions_example, pathwayType = "KEGG",
    max_homology_bitscore = 60, precision_window = 400, exclude_pathways = "blacklist")

example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )



source("https://bioconductor.org/biocLite.R")
biocLite("STRINGdb")

dataDecoratedByString <- decorateByStringData(
    chebiMoleculesDecoratedByReactomeDf = decReac[decReac$root %in% c("CHEBI:15756","CHEBI:35465"),]
)
cat("\014")
diff_exp_example2 <- data.frame("gene" = dfElement[,"genesSymbolsFromEnsemble"][[1]])

decReac[decReac$root %in% c("CHEBI:15756","CHEBI:35465"),][1,"ensembleIds"][[1]]


diff_exp_example2 <- data.frame("translate" = decReac[,"ensembleIds"][[1]])
example2_mapped <- string_db$map( diff_exp_example2, "translate", removeUnmappedRows = TRUE )
example2_mapped[,"STRING_id"]
res <- plyr::aaply(.data = example2_mapped[,"STRING_id"],
            .margins = 1, .fun = function(x) {
                z <- igraph::neighbors(ggg, x)
                print(length(z$name))
                length(z$name)
            })


plyr::ddply(.data = example2_mapped, .variables = c("STRING_id"), .fun = function(r) {
    print(r)
    r
})

ggg <- string_db$get_graph()
ssssaaaa<-plyr::ddply(.data = example2_mapped, .variables = c("STRING_id"), .fun = function(r) {
    data.frame("res" = I(list(igraph::neighbors(ggg, r[,"STRING_id"])$name)))
})

unique(unlist(ssssaaaa[,"res"]))
intersect(ssssaaaa[,"res"][[1]],ssssaaaa[,"res"][[2]])


cat("\014")
unique(unlist(ssssaaaa[,"res"]))

interSect <- unique(ssssaaaa[,"res"][[1]])
for(i in length(ssssaaaa[,"res"])) {
    interSect <- intersect(interSect, ssssaaaa[,"res"][[i]])
}
ssssaaaa[,"res"][[2]]


totalVector <- c()
for() {

}

str(ssssaaaa)

length(example2_mapped[,"STRING_id"])

decReac[decReac$root %in% c("CHEBI:15756","CHEBI:35465"),][1,"ensembleIds"][[1]]

chebiIdsToReactomePathways = decReac[decReac$root %in% c("CHEBI:15756","CHEBI:35465"),]
stringOrganismId = 9606
stringDbVersion = "10"
idsColumnName = 'ontologyId'
rootColumnName = 'root'
listOfEnsembleIdColumnName = 'ensembleIds'
    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(idsColumnName, rootColumnName)
    }
    string_db <- STRINGdb$new(
        version = stringDbVersion,
        species = stringOrganismId,
        input_directory = path.expand("~"))
    chebiIdsToRealReactomePathways <- chebiIdsToReactomePathways[!chebiIdsToReactomePathways[idsColumnName] == "", ]
    # dfWithString <- plyr::ddply(.data = chebiIdsToRealReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
    dfElement <- chebiIdsToRealReactomePathways
        print(dfElement)
        returnNeighbourVector <- character(length = 0)
        stringGenesSymbols <- character(length = 0)
        if (0 == length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
        } else {
            proteinIds <- getEnsemblProteinsIdsBaseOnEnsemblGensIdsUsingMyGenePackage(
                dfElement[1,listOfEnsembleIdColumnName], organismTaxonomyId = stringOrganismId
            )
            stringId1 <- string_db$mp(proteinIds)
            returnNeighbourVector <- string_db$get_neighbors(stringId1)
            string_db$get_neighbors()
            ensembleIdsFromStringDb <- mapFromStringIdsToEnsembleIds(returnNeighbourVector)
            stringGenesSymbols <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDb, organismTaxonomyId = stringOrganismId)
        }
        dffff <- data.frame('ensembleIds' = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(returnNeighbourVector))),
                            'stringGenesSymbols' = I(list(unique(stringGenesSymbols))) )
        dffff

    # })
    dfWithString




























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


