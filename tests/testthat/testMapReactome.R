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

test_that("NewOnionApiWorkflow test", {
    #given
    print('*********given*********')
    clasteredSmallMolecules <- ONION::clasterUsingOntology(pathToFile = "/home/koralgooll/doktorat/Rpackages/ONION/example/nm-lipidomics-valid.txt",
                                header = TRUE, ontologyRepresentatnion = ONION::createFirstExistsInReactomeChebiOntology)
    mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(clasteredSmallMolecules)
    chebiIdsToReactomePathways <- ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules, organismTaxonomyId = '9606')
    chebiIdsToReactomePathways$ensembleIds[[11]]
    chebiIdsToReactomePathwaysAndToStringNeighbours <- ONION::getStringNeighbours(chebiIdsToReactomePathways)
})

test_that("shouldGoThroughONIONAPIWorkflow test", {
    #given
    print('*********given*********')
    smallMolecules <- clusterSmallMolecules("../../example/smallMolecules.txt")
    smallMolecules <- clusterSmallMolecules("./example/smallMolecules.txt")
    mergeSM <- mergeChEBIOntologyWithChildFavoring(smallMolecules)
    #ID mapping. Reactome <-> TaxonId.
    ms <- mapReactomePathways(mergeSM)
    mp <- getStringNeighbours(ms)

    #when
    msPseudoC <- showPseudoClusteringResultsOnGens(ms)
    mpPseudoC <- showPseudoClusteringResults(mp)


    resS <- makeCCAOnData("15756", msPseudoC$`15756`$hgnc_symbol, "D:/doktorat/repositories/ONION/example/nm-transcriptomics.txt",
                         "D:/doktorat/repositories/ONION/example/nm-lipidomics.txt")

    #Matrix created from ONION and nm-transcriptomics.txt is too big for CCA. This will throw error.
    resP <- makeCCAOnData("15756", mpPseudoC$`15756`$hgnc_symbol, "D:/doktorat/repositories/ONION/example/nm-transcriptomics.txt",
                         "D:/doktorat/repositories/ONION/example/nm-lipidomics.txt")
    #Smaller set givs results. :)
    resP <- makeCCAOnData("15756", mpPseudoC$`15756`$hgnc_symbol, "D:/doktorat/repositories/ONION/example/transTest.txt",
                          "D:/doktorat/repositories/ONION/example/nm-lipidomics.txt")

    #then
    expect_that( resS, is_a("list") )
})

