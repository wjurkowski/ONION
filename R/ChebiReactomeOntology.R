
firstExistsInReactomeChebiOntology <- function(baseData, rootColumnName) {
    # print("firstExistsInReactomeChebiOntology")
    # print(as.character(baseData[, rootColumnName]))
    firstExistsInReactomeChebiOntologyDataFrame <- ldply(as.character(baseData[, rootColumnName]), function(dataRow) {
        # print("firstExistsInReactomeChebiOntology:BEGINING")
        child <- 0
        parent <- 0
        # str(dataRow)
        # print(dataRow)
        # print()
        if (checkIfPathwayIdExistsForChEBIId(strsplit(as.character(dataRow), ":")[[1]][2])) {
            child <- 1
            parent <- 1
        } else {
            # print("firstExistsInReactomeChebiOntology:CHILDREN-PARENTS")
            children <- getChEBIOntologyChildren(strsplit(as.character(dataRow), ":")[[1]][2])
            parents <- getChEBIOntologyParents(strsplit(as.character(dataRow), ":")[[1]][2])
            # print(children)
            # print(parents)
            if (nrow(children) >= 1) {
                for (i in 1:nrow(children)) {
                    # print(children[i, rootColumnName])
                    # print(as.character(children[i, rootColumnName]))
                    if (checkIfPathwayIdExistsForChEBIId(as.character(children[i, ]$chebiId))) {
                        child <- paste("CHEBI", as.character(children[i, ]$chebiId), collapse = NULL, sep = ":")
                        break
                    }
                }
            }
            if (nrow(parents) >= 1) {
                for (i in 1:nrow(parents)) {
                    # print(parents[i, rootColumnName])
                    # print(as.character(parents[i, rootColumnName]))
                    if (checkIfPathwayIdExistsForChEBIId(as.character(parents[i, ]$chebiId))) {
                        parent <- paste("CHEBI", as.character(parents[i, ]$chebiId), collapse = NULL, sep = ":")
                        break
                    }
                }
            }
        }
        data.frame('child' = as.character(child),
                   'root' = as.character(dataRow),
                   'parent' = as.character(parent))
    })
    firstExistsInReactomeChebiOntologyDataFrame
}

# NEW PUBLIC API:
mergeChEBIOntologyWithChildFavoring <- function(chebiOntologyDataFrame, rootColumnName = 'root') {
    mergeOntologyDataFrame <- ddply(chebiOntologyDataFrame, c(rootColumnName), function(dataFrameRow){
        ontologyId <- ""
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
        data.frame('root' = dataFrameRow$root, 'ontologyId' = ontologyId, 'whoWins' = whoWins)
    })
    mergeOntologyDataFrame
}