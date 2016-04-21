
cat('\014')
chebiIdsToReactomePathwaysAndToStringNeighbours


makeCanonicalCorrelationAnalysis  <- function(vectorOfChebiIds, vectorOfHgncSymbols,
                          pathToFileWithTranscriptomicsData, pathToFileWithLipidomicsData) {
    transcriptomicsData <- readWithoutDuplicates(pathToFileWithTranscriptomicsData)
    lipidomicsData <- readWithoutDuplicates(pathToFileWithLipidomicsData)
    X <- as.matrix(t(transcriptomicsData))
    Y <- as.matrix(t(lipidomicsData))
    matchedGensData <- match(vectorOfHgncSymbols, colnames(X))
    matchedChebiData <- match(vectorOfChebiIds, colnames(Y))
    factorOfMatchedGensData <- factor(matchedGensData)
    factorOfMatchedChebiData <- factor(matchedChebiData)
    machedX <- X[,as.numeric(levels(factorOfMatchedGensData))]
    machedY <- Y[,as.numeric(levels(factorOfMatchedChebiData))]

    if (is.numeric(machedX) && is.numeric(machedX)) {
        if (is.matrix(machedX)) {
        } else {
            machedX <- matrix(machedX)
        }
        if (is.matrix(machedY)) {
        } else {
            machedY <- matrix(machedY)
        }
    }

    ccaFromyacca <- tryCatch(
        {
            cca(machedX, machedY)
        },
        error=function(cond) {
            message("ONION - Included CCA can not solve task.")
            message("Original message (yacca):")
            message(cond)
            # Choose a return value in case of error
            return(NA)
        },
        warning=function(cond) {
            message("ONION - Included CCA present warning.")
            message("Original message (yacca):")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
            message("ONION - CCA (yacca) finished with success.")
        }
    )

    #doubleCCA = list(CCA = ccaFromCCA, yacca = ccaFromyacca)
    #doubleCCA
    ccaFromyacca
}

