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

test_that("shouldGoThroughONIONAPIWorkflow test", {
    #given
    print('*********given*********')
    smallMolecules <- clusterSmallMolecules("../../example/smallMolecules.txt")
    mergeSM <- mergeChEBIOntologyWithChildFavoring(smallMolecules)
    #ID mapping. Reactome <-> TaxonId.
    ms <- mapReactomePathways(mergeSM)
    mp <- getStringNeighbours(ms)

    #when
    pseudoClustering <- showPseudoClusteringResults(mp)

    #then
    expect_that( pseudoClustering, is_a("list") )
    #pseudoClustering[["36023"]]
})

