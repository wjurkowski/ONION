#Set working directory. Solve relative path problem in tests.
#TODO Add different behavior on CHECK and TEST.
setwd(paste(getwd(), "/.."))
print(getwd())

source("../R/mapReactome.R", chdir = TRUE)


set_up <- function() {
    source("R/mapReactome.R")
}

tear_down <- function() {
    cat("\014")
}


test_that("clusterSmallMolecules test", {
    #given

    #when
    mr <- clusterSmallMolecules("C:/HOME/ONIONpackage/ONION/R/smallMolecules.txt")
    mr
    #then
    expect_that( mr, is_a("data.frame") )
})

test_that("mapReactomePathways test", {
    #given
    smallMolecules <- clusterSmallMolecules("C:/HOME/ONIONpackage/ONION/R/smallMolecules.txt")

    #when
    margeSM <- margeChEBIOntologyWithChildFavoring(smallMolecules)
    margeSM
    ms <- mapReactomePathways(margeSM, "HSA")
    ms

    #then
    expect_that( ms, is_a("list") )
})

test_that("getStringNeighbours test", {
    #given
    smallMolecules <- clusterSmallMolecules("C:/HOME/ONIONpackage/ONION/R/smallMolecules.txt")
    margeSM <- margeChEBIOntologyWithChildFavoring(smallMolecules)
    ms <- mapReactomePathways(margeSM, "HSA")

    #when
    mp <- getStringNeighbours(ms)
    mp

    #then
    expect_that( mp, is_a("list") )
})

test_that("showPseudoClusteringResults test", {
    #given
    smallMolecules <- clusterSmallMolecules("C:/HOME/ONIONpackage/ONION/R/smallMolecules.txt")
    margeSM <- margeChEBIOntologyWithChildFavoring(smallMolecules)
    #ID mapping. Reactome <-> TaxonId.
    ms <- mapReactomePathways(margeSM, "HSA")
    mp <- getStringNeighbours(ms)

    #when
    pseudoClustering <- showPseudoClusteringResults(mp)

    #then
    expect_that( pseudoClustering, is_a("list") )
    pseudoClustering[["36023"]]
})

