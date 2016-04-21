# TRUE help
readWithoutDuplicates <- function(filePath, header) {
    transcriptomics <- read.table(filePath, header = header, row.names = NULL);
    transcriptomicsWithoutDuplicates <- data.frame(transcriptomics[!duplicated(transcriptomics[1]), ], row.names = 1)
}


#TODO move to utils
#priver util methods
removeEmptyLists <- function(dirtyList) {
    recurention <- FALSE
    if (length(dirtyList) == 0) {
    } else {
        for (i in 1:length(dirtyList)) {
            if (length(dirtyList[[i]]) == 0) {
                dirtyList[[i]] <- NULL
                recurention <- TRUE
                break
            }
        }
    }
    if (recurention) {
        dirtyList <- removeEmptyLists(dirtyList)
    }
    dirtyList
}

flattenList <- function(listToFlatten) {
    lapply(listToFlatten, function(listElement){
        if (length(listElement) != 0) {
            listToReturn <- list()
            for (i in 1:length(listElement)) {
                if (length(listElement[[i]]) != 0) {
                    for (j in 1:length(listElement[[i]])) {
                        listToReturn <- c(listToReturn, listElement[[i]][[j]])
                    }
                }
            }
            #as.vector(listToReturn)
            as.character(listToReturn)
            #listToReturn
        }
    })
}

removeDuplicatesInListVectors <- function(listOfVectors) {
    listOfVectorsWithoutDuplicates <- lapply(listOfVectors, function(listElement){
        unique(x = listElement)
    })
    listOfVectorsWithoutDuplicates
}

addNames <- function(ChEBIList, flattenedList) {
    for (i in 1:length(flattenedList)) {
        names(flattenedList)[i] <- ChEBIList[i,]$ChEBI
    }
    flattenedList
}

removeEmptyElementsOnListOfCharsVectors <- function(listOfVectors) {
    listOfVectorsWithoutEmptyChars <- lapply(listOfVectors, function(listElement){
        clearEnsemblePeptideIds <- listElement[listElement != ""]
        clearEnsemblePeptideIds
    })
    listOfVectorsWithoutEmptyChars
}

checkIfValueCanBeProcessed <- function(value) {
    if (isLengthZero(value) || isCharNA(value)) {
        FALSE
    } else {
        TRUE
    }
}

isCharNA <- function(value) {
    if ("NA" == value) {
        TRUE
    } else {
        FALSE
    }
}

isLengthZero <- function(value) {
    if (length(value) == 0) {
        TRUE
    } else {
        FALSE
    }
}