
firstExistsInReactomeChebiOntology <- function(baseData) {
    firstExistsInReactomeChebiOntologyDataFrame <- ldply(baseData$ChEBI, function(dataRow) {
        child <- 0
        parent <- 0
        if (checkIfPathwayIdExistsForChEBIId(dataRow)) {
            child <- 1
            parent <- 1
        } else {
            children <- getChEBIOntologyChildren(dataRow)
            parents <- getChEBIOntologyParents(dataRow)
            for (i in 1:nrow(children)) {
                if (checkIfPathwayIdExistsForChEBIId(as.character(children[i, ]$chebiId))) {
                    child <- as.character(children[i, ]$chebiId)
                    break
                }
            }
            for (i in 1:nrow(parents)) {
                if (checkIfPathwayIdExistsForChEBIId(as.character(parents[i, ]$chebiId))) {
                    parent <- as.character(parents[i, ]$chebiId)
                    break
                }
            }
        }
        data.frame('child' = as.numeric(child),
                   'root' = as.numeric(dataRow),
                   'parent' = as.numeric(parent))
    })
    firstExistsInReactomeChebiOntologyDataFrame
}

# NEW PUBLIC API:
mergeChEBIOntologyWithChildFavoring <- function(chebiOntologyDataFrame) {
    mergeOntologyDataFrame <- ddply(chebiOntologyDataFrame, .(root), function(dataFrameRow){
        ontologyId <- 0
        whoWins <- 'N'
        if (1 == dataFrameRow$child && 1 == dataFrameRow$parent) {
            ontologyId <- dataFrameRow$root
            whoWins <- 'R'
        } else if (0 != dataFrameRow$child) {
            ontologyId <- dataFrameRow$child
            whoWins <- 'C'
        } else if (0 != dataFrameRow$parent) {
            ontologyId <- dataFrameRow$parent
            whoWins <- 'P'
        }
        data.frame('ontologyId' = ontologyId, 'root' = dataFrameRow$root, 'whoWins' = whoWins)
    })
    mergeOntologyDataFrame
}