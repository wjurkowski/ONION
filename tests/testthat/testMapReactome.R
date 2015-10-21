#Set working directory. Solve relative path problem in tests.
#TODO Add different behavior on CHECK and TEST.
setwd(paste(getwd(), "/.."))
print(getwd())

source("../R/mapReactome.R", chdir = TRUE)

test_that("clusterSmallMolecules test", {
    #given

    #whene
    mr <- clusterSmallMolecules("C:/HOME/ONIONpackage/ONION/R/smallMolecules.txt")()

    #then
    expect_that( mr, is_a("data.frame") )
})
