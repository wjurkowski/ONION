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
    smallMolecules <- clusterSmallMolecules("D:/doktorat/repositories/ONION/example/smallMolecules.txt")
    margeSM <- margeChEBIOntologyWithChildFavoring(smallMolecules)
    #ID mapping. Reactome <-> TaxonId.
    ms <- mapReactomePathways(margeSM)
    mp <- getStringNeighbours(ms)

    #when
    pseudoClustering <- showPseudoClusteringResults(mp)

    #then
    expect_that( pseudoClustering, is_a("list") )
    #pseudoClustering[["36023"]]
})

