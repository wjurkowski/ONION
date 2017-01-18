
biopaxExporter <- function(level = "level3", eventId) {
    url <- paste("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter", level, eventId, sep = "/")
    response <- getForm(url, .opts = list(httpheader = c('Content-Type' = 'application/json')))
    bioPAXList <- xmlToList(response[[1]])
}

frontPageItems <- function(speciesName="homo+sapiens") {
    url <- paste("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/frontPageItems", speciesName, sep = "/")
    response <- getForm(url, .opts = list(httpheader = c('Content-Type' = 'application/json')))
    frontPageItemsList <- fromJSON(response[[1]])
}

highlightPathwayDiagram <- function(pathwayId) {
    #TODO Implement.
}

listByQuery <- function(className="", query, ...) {
    #TODO Implement.
}

pathwayDiagram <- function(pathwayId, returnedFormat="PNG") {
    #TODO Implement.
}

pathwayHierarchy <- function(speciesName="homo+sapiens") {
    #TODO Implement.
}

pathwayParticipants <- function(pathwayId) {
    #TODO Implement.
}

pathwayComplexes <- function(pathwayId) {
    #TODO Implement.
}

queryById <- function(className="Pathway", dbId) {
    url <- paste("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/queryById", className, dbId, sep = "/")
    response <- getForm(url, .opts = list(httpheader = c('Content-Type' = 'application/json')))
    dbObjectList <- fromJSON(response[[1]])
}

queryByIds <- function(className="Pathway", dbId) {
    url <- paste("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/queryByIds", className, sep = "/")
    response <- postForm(url, id=dbId, style = "POST", .opts = list(httpheader = c('Content-Type' = 'text/plain')))
    dbObjectsList <- fromJSON(response[[1]])
}

queryHitPathways <- function() {
    #TODO It is more sophisticated case.
    #body="PPP2R1A,CEP192,AKAP9,CENPJ,CEP290,DYNC1H1"
}

#TODO It works only with one ID, but it should work with few coma separated.
pathwaysForEntitie <- function(physicalEntityDatabaseId) {
    url <- paste("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/pathwaysForEntities")
    response <- postForm(url, "ID"=physicalEntityDatabaseId, style = "POST", .opts = list(httpheader = c('Content-Type' = 'application/json')))
    pathwaysList <- fromJSON(response[[1]])
}

queryReviewedPathways <- function(personId) {
    #TODO Implement.
}

queryPeopleByName <- function(name) {
    #TODO Implement.
}

queryPeopleByEmail <- function(email) {
    #TODO Implement.
}

speciesList <- function() {
    response <- getForm("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/speciesList",
                        .opts = list(httpheader = c('Content-Type' = 'application/json')))
    speciesList <- fromJSON(response[[1]])
}

sbmlExporter <- function(eventId) {
    #TODO Implement.
}

getReferenceMolecules <- function() {
    response <- getForm("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/getReferenceMolecules",
                        .opts = list(httpheader = c('Content-Type' = 'application/json')))
    referenceMolecules <- stringResponseBodyToDataFrame(response[[1]], "\t", "reactomeId", "chebiId")
}

getDiseases <- function() {
    response <- getForm("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/getDiseases",
                        .opts = list(httpheader = c('Content-Type' = 'application/json')))
    diseases <- stringResponseBodyToDataFrame(response[[1]], "\t", "rectomeId", "diseaseOntologyId")
}

getUniProtRefSeqs <- function() {
    response <- getForm("http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/getUniProtRefSeqs",
                        .opts = list(httpheader = c('Content-Type' = 'application/json')))
    uniProtRefSeqs <- stringResponseBodyToDataFrame(response[[1]], "\t", "reactomeId", "uniprotId")
}

#GlobalForPerformance
.onLoad <- function(libname, pkgname) {
    chEBIToReactomeLowestLevel <<- read.table(paste(find.package("ONION"), "resources/ChEBI2Reactome.txt", sep = "/"))
    Ensembl2ReactomeLowestLevel <<- read.table(paste(find.package("ONION"), "resources/Ensembl2Reactome.txt", sep = "/"))
    UniProt2ReactomeLowestLevel <<- read.table(paste(find.package("ONION"), "resources/UniProt2Reactome.txt", sep = "/"))
}

#Additional (local) API methods
checkIfPathwayIdExistsForChEBIId <- function(ChEBIId) {
    print("checkIfPathwayIdExistsForChEBIId")
    print(ChEBIId)
    if (0 == length(ChEBIId)) {
        FALSE
    } else {
        matchingChEBIToReactome <- chEBIToReactomeLowestLevel[grep(ChEBIId, chEBIToReactomeLowestLevel$V1),]
        if (nrow(matchingChEBIToReactome)) {
            TRUE
        } else {
            FALSE
        }
    }

}

# OLD
getEnsemblIdsForPathway <- function(pathwayId) {
    matchingEnsemblsToPathway <- Ensembl2ReactomeLowestLevel[grep(pathwayId, Ensembl2ReactomeLowestLevel$V2),]
    as.list(as.character(matchingEnsemblsToPathway$V1))
}

# NEW
getEnsemblIdsForPathwayId <- function(pathwayId) {
    matchingEnsemblsToPathway <- Ensembl2ReactomeLowestLevel[grep(pathwayId, Ensembl2ReactomeLowestLevel$V2),]
    as.character(matchingEnsemblsToPathway$V1)
}

# NEW
getUniProtIdsForPathwayId <- function(pathwayId) {
    matchingEnsemblsToPathway <- UniProt2ReactomeLowestLevel[grep(pathwayId, UniProt2ReactomeLowestLevel$V2),]
    as.character(matchingEnsemblsToPathway$V1)
}

# TOTALY NEW
getEnsemblIdsForPathwayIds <- function(vectorOfPAthwayIds) {
    ensembleIdsList <- alply(.data = vectorOfPAthwayIds, .fun = function(vectorElement) {
        getEnsemblIdsForPathwayId(vectorElement)
    }, .margins = 1)
    ensembreIdsVector <- as.character(unique(unlist(ensembleIdsList, recursive = TRUE)))
    ensembreIdsVector
}

# TOTALY NEW
getUniProtIdsForPathwayIds <- function(vectorOfPAthwayIds) {
    uniProtIdsList <- alply(.data = vectorOfPAthwayIds, .fun = function(vectorElement) {
        getUniProtIdsForPathwayId(vectorElement)
    }, .margins = 1)
    uniProtIdsVector <- as.character(unique(unlist(uniProtIdsList, recursive = TRUE)))
    uniProtIdsVector
}

#OLD
getPathwaysIdsForChEBI <- function(ChEBIId) {
    matchingChEBIToReactome <- chEBIToReactomeLowestLevel[grep(ChEBIId, chEBIToReactomeLowestLevel$V1),]
    as.list(as.character(matchingChEBIToReactome$V2))
}

#NEW
getPathwaysIdsForChebiId <- function(ChEBIId) {
    # TODO: Validation on ChEBIId!
    if (0 == ChEBIId) {
        returnValue <- character(length = 0)
    } else {
        matchingChEBIToReactome <- chEBIToReactomeLowestLevel[grep(ChEBIId, chEBIToReactomeLowestLevel$V1),]
        returnValue <- as.character(matchingChEBIToReactome$V2)
    }
    returnValue
}

# OLD
getPathwaysIdsForChEBIAndOrganismCode <- function(ChEBIId, organismCode) {
    dirtyListOfPathwaysIds <- getPathwaysIdsForChebiId(ChEBIId)
    clearListOfPathwaysIds <- lapply(dirtyListOfPathwaysIds, function(listElement) {
        if (strsplit(listElement, "-")[[1]][2] == organismCode) {
            #strsplit(listElement, "-")[[1]][3]
            listElement
        } else {
            NA
        }

    })
    clearListOfPathwaysIds
    withoutDuplicates <- removeDuplicates(clearListOfPathwaysIds)
    removeNA(withoutDuplicates)
}

# NEW
getPathwaysIdsForChebiUnderOrganism <- function(ChEBIId, organismCode) {
    dirtyVectorOfPathwaysIds <- getPathwaysIdsForChebiId(ChEBIId)
    vectorsOfPathwayIdsUnderOrganism <- aaply(.data = dirtyVectorOfPathwaysIds, .fun = function(vectorElement) {
        if (strsplit(vectorElement, "-")[[1]][2] == organismCode) {
            vectorElement
        } else {
            NA
        }
    }, .margins = 1)
    cleanVectorOfPathwaysIds <- as.character(unique(vectorsOfPathwayIdsUnderOrganism[!is.na(vectorsOfPathwayIdsUnderOrganism)]))
    cleanVectorOfPathwaysIds
}

removeDuplicates <- function(dirtyList) {
    duplicationsInDirtyList <- as.list(duplicated(dirtyList, incomparables = FALSE))
    noDuplicationVector <- list()
    for(i in 1:length(duplicationsInDirtyList)) {
        if (!duplicationsInDirtyList[[i]]) {
            noDuplicationVector <- c(noDuplicationVector, dirtyList[i])
        }
    }
    noDuplicationVector
}

removeNA <- function(dirtyList) {
    recurention <- FALSE
    for (i in 1:length(dirtyList)) {
        if (is.na(dirtyList[[i]])) {
            dirtyList[[i]] <- NULL
            recurention <- TRUE
            break
        }
    }
    if (recurention) removeNA(dirtyList)
    dirtyList
}

getUsablePathwaysIdsForChEBIForAllOrganisms <- function(ChEBIId) {
    dirtyListOfPathwaysIds <- getPathwaysIdsForChEBI(ChEBIId)
    clearListOfPathwaysIds <- lapply(dirtyListOfPathwaysIds, function(listElement) {
        strsplit(listElement, "-")[[1]][3]
    })
    clearListOfPathwaysIds
    removeDuplicates(clearListOfPathwaysIds)
}

#private methods
stringResponseBodyToDataFrame <- function(stringResponseBody, sep = "\t", ...) {
    responseConnection <- textConnection(stringResponseBody)
    dataFrame <- read.table(responseConnection, sep=sep, col.names=c(...),
                            fill=FALSE, strip.white=TRUE)
}
