
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