testSuite <- function() {
    # Before tests
    pathToFileWithChebiIds <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")
    baseData <- read.table(pathToFileWithChebiIds, header = TRUE)

    testthat::test_that("one row DF clustering", {
        # Given
        firstRow <- baseData[1,]

        # When
        clusteredSmallMolecule <- ONION::clusterUsingOntology(
            chebiIdsDataFrame = firstRow,
            rootColumnName = "ChEBI",
            ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology
        )

        # Then
        testthat::expect_identical(
            clusteredSmallMolecule,
            data.frame(child = as.factor(c(1)), root = c("CHEBI:28875"), parent = as.factor(c(1)))
        )

        # Uncomment this please, it doesn't pass test.
        # testthat::expect_identical(
        #     clusteredSmallMolecule,
        #     data.frame(child = c(1), root = c("CHEBI:28875"), parent = as.factor(c(1)))
        # )
    })

    testthat::test_that("find child while clustering", {
        # Given
        rowWithChild <- baseData[5,]

        # When
        clusteredSmallMolecule <- ONION::clusterUsingOntology(
            chebiIdsDataFrame = rowWithChild,
            rootColumnName = "ChEBI",
            ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology
        )

        # Then
        testthat::expect_identical(
            clusteredSmallMolecule$child,
            as.factor(c("CHEBI:15541"))
        )

        # Uncomment this please, it doesn't pass test.
        # testthat::expect_identical(
        #     clusteredSmallMolecule$child,
        #     as.factor(c("CHEBI:11111"))
        # )
    })

    testthat::test_that("find child while clustering", {
        # Given
        NotDataFrame <- character(0)

        # When and Then
        testthat::expect_error(ONION::clusterUsingOntology(
            chebiIdsDataFrame = NotDataFrame,
            rootColumnName = "ChEBI",
            ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology
        ), "incorrect number of dimensions")

        # Uncomment this please, it doesn't pass test.
        # testthat::expect_error(ONION::clusterUsingOntology(
        #     chebiIdsDataFrame = NotDataFrame,
        #     rootColumnName = "ChEBI",
        #     ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology
        # ), "not real error message")
    })

}

testSuite()
