
createFirstExistsInReactomeChebiOntology <- function(baseData) {
    firstExistsInReactomeChebiOntologyDataFrame <- ldply(baseData$ChEBI, function(dataRow) {
        child <- 0
        parent <- 0
        if (ReactomeAPI::checkIfPathwayIdExistsForChEBIId(dataRow)) {
            child <- 1
            parent <- 1
        } else {
            children <- ChebiAPI::getChEBIOntologyChildren(dataRow)
            parents <- ChebiAPI::getChEBIOntologyParents(dataRow)
            for (i in 1:nrow(children)) {
                if (ReactomeAPI::checkIfPathwayIdExistsForChEBIId(as.character(children[i, ]$chebiId))) {
                    child <- as.character(children[i, ]$chebiId)
                    break
                }
            }
            for (i in 1:nrow(parents)) {
                if (ReactomeAPI::checkIfPathwayIdExistsForChEBIId(as.character(parents[i, ]$chebiId))) {
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
        if (1 == dataFrameRow$child && 1 == dataFrameRow$parent) {
            ontologyId <- dataFrameRow$root
        } else if (0 != dataFrameRow$child) {
            ontologyId <- dataFrameRow$child
        } else if (0 != dataFrameRow$parent) {
            ontologyId <- dataFrameRow$parent
        }
        data.frame('ontologyId' = ontologyId)
    })
    mergeOntologyDataFrame['ontologyId']
}