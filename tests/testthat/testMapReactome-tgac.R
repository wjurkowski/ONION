#Set working directory. Solve relative path problem in tests.
#TODO Add different behavior on CHECK and TEST.
#setwd(paste(getwd(), "/.."))
#print(getwd())

#setup
install.packages("/Users/jurkowsw/apps/ONION/builds/sourcePackages/ONION_0.0.1.tar.gz", repos = NULL, type = 'source')
library(ONION)

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
    #smallMolecules <- clusterSmallMolecules("./example/smallMolecules.txt")
    mergeSM <- mergeChEBIOntologyWithChildFavoring(smallMolecules)
    #ID mapping. Reactome <-> TaxonId.
    ms <- mapReactomePathways(mergeSM)
    mp <- getStringNeighbours(ms)

    #when
    msPseudoC <- showPseudoClusteringResultsOnGens(ms)
    mpPseudoC <- showPseudoClusteringResults(mp)


    resS <- makeCCAOnData("15756", msPseudoC$`15756`$hgnc_symbol, "../../example/nm-transcriptomics.txt",
                         "../../example/nm-lipidomics.txt")

    #Matrix created from ONION and nm-transcriptomics.txt is too big for CCA. This will throw error.
    resP <- makeCCAOnData("15756", mpPseudoC$`15756`$hgnc_symbol, "../../example/nm-transcriptomics.txt",
                         "../../example/nm-lipidomics.txt")
  
    
    #Smaller set givs results. :)
    resP <- makeCCAOnData("15756", mpPseudoC$`15756`$hgnc_symbol, "D:/doktorat/repositories/ONION/example/transTest.txt",
                          "D:/doktorat/repositories/ONION/example/nm-lipidomics.txt")

    #then
    expect_that( resS, is_a("list") )
})

