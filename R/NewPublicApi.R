
#GlobalForPerformance
setUpReactomeMapping <- function(ChEBI2ReactomeFileURL, Ensembl2ReactomeFileURL, UniProt2ReactomeFileURL) {

    ChEBI2ReactomeFileURL <- RCurl::getURL(ChEBI2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    chEBIToReactomeLowestLevel <<- read.table(
        textConnection(ChEBI2ReactomeFileURL),
        header = FALSE
    )

    Ensembl2ReactomeFileURL <- RCurl::getURL(Ensembl2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    Ensembl2ReactomeLowestLevel <<- read.table(
        textConnection(Ensembl2ReactomeFileURL),
        header = FALSE
    )

    UniProt2ReactomeFileURL <- RCurl::getURL(UniProt2ReactomeFileURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    UniProt2ReactomeLowestLevel <<- read.table(
        textConnection(UniProt2ReactomeFileURL),
        header = FALSE
    )
}

decorateByReactomeData <- function(chebiMoleculesDf, chebiIdsColumnName, organismTaxonomyId = '9606') {

    clusteredSmallMolecules <- OmicsON::clusterUsingOntology(
        chebiIdsDataFrame = chebiMoleculesDf,
        rootColumnName = chebiIdsColumnName,
        ontologyRepresentatnion = OmicsON::firstExistsInReactomeChebiOntology)

    mergedSmallMolecules <- OmicsON::mergeChEBIOntologyWithChildFavoring(
        clusteredSmallMolecules,
        rootColumnName = 'root')

    OmicsON::mapReactomePathwaysUnderOrganism(
        chebiOntologyIds = mergedSmallMolecules[, c("ontologyId", "root"), drop = FALSE],
        organismTaxonomyId = organismTaxonomyId,
        idsColumnName = "ontologyId",
        rootColumnName = "root")
}

decorateByStringDbData <- function(chebiIdsToReactomePathways, listOfEnsembleIdColumnName,
                                   stringOrganismId = '9606', stringDbVersion = "10",
                                   idsColumnName = 'ontologyId', rootColumnName = 'root') {
    stringOrganismId <- as.numeric(stringOrganismId)

    if (is.null(rootColumnName)) {
        columnsUseInIteration <- c(idsColumnName)
    } else {
        columnsUseInIteration <- c(rootColumnName, idsColumnName)
    }

    string_db <- STRINGdb$new(
        version = stringDbVersion,
        species = stringOrganismId,
        input_directory = path.expand("~"))
    dfWithString <- ddply(.data = chebiIdsToReactomePathways, columnsUseInIteration, .fun = function(dfElement) {
        extendedByStringAsVector <- character(length = 0)
        stringGenesSymbolsExpand <- character(length = 0)
        stringGenesSymbolsNarrow <- character(length = 0)
        if (0 != length(dfElement[1, listOfEnsembleIdColumnName][[1]])) {
            toTranslate <- data.frame("translate" = dfElement[1,listOfEnsembleIdColumnName][[1]])
            translated <- string_db$map( toTranslate, "translate", removeUnmappedRows = TRUE )
            stringGraph <- string_db$get_graph()

            extendedByString <- plyr::ddply(.data = translated, .variables = c("STRING_id"), .fun = function(r) {
                reksultDataFrame <- data.frame("res" = I(list(c(""))))
                try((function() {
                    reksultDataFrame <<- data.frame("res" = I(list(igraph::neighbors(stringGraph, r[,"STRING_id"])$name)))
                })())
                reksultDataFrame
            })

            extendedByStringAsVector <- unique(unlist(extendedByString[,"res"]))
            ensembleIdsFromStringDbExpand <- mapFromStringIdsToEnsembleIds(extendedByStringAsVector)
            stringGenesSymbolsExpand <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDbExpand,
                                                                                           organismTaxonomyId = stringOrganismId)

            interSect <- character(length = 0)
            if (1 >= length(extendedByString[,"res"])) {
                interSect <- character(length = 0)
            } else {
                interSect <- unique(extendedByString[,"res"][[1]])
                for(i in length(extendedByString[,"res"])) {
                    interSect <- intersect(interSect, extendedByString[,"res"][[i]])
                }
            }
            ensembleIdsFromStringDbNarrow <- mapFromStringIdsToEnsembleIds(interSect)
            stringGenesSymbolsNarrow <- getSymbolsBaseOnEnsemblPeptidIdsUsingMyGenePackage(ensembleIdsFromStringDbNarrow,
                                                                                           organismTaxonomyId = stringOrganismId)
        }

        dffff <- data.frame(listOfEnsembleIdColumnName = dfElement[1,listOfEnsembleIdColumnName][1],
                            'stringIds' = I(list(unique(extendedByStringAsVector))),
                            'stringGenesSymbolsExpand' = I(list(unique(stringGenesSymbolsExpand))),
                            'stringGenesSymbolsNarrow' = I(list(unique(stringGenesSymbolsNarrow))))
        dffff
    })
    names(dfWithString)[3] <- listOfEnsembleIdColumnName
    dfWithString
}
