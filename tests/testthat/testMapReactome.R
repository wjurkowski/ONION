#Set working directory. Solve relative path problem in tests.
#TODO Add different behavior on CHECK and TEST.
setwd(paste(getwd(), "/.."))
print(getwd())

source("../R/mapReactome.R", chdir = TRUE)

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
    smallMolecules
    margeSM <- margeChEBIOntologyWithChildFavoring(smallMolecules)
    margeSM
    ms <- mapReactomePathways(margeSM, "HSA")
    ms
    ms[["28364"]]

    #then
    expect_that( mr, is_a("list") )
})
