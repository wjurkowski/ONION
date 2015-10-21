library(RCurl)
library(RJSONIO)
library(XML)

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
chEBIToReactomeLowestLevel <- read.table("resources/ChEBI2Reactome.txt")
#Additional (local) API methods
checkIfPathwayIdExistsForChEBIId <- function(ChEBIId) {
    #chEBIToReactomeLowestLevel <- read.table("resources/ChEBI2Reactome.txt")
    #chEBIToReactomeLowestLevel[chEBIToReactomeLowestLevel$V1=="28842",]
    matchingChEBIToReactome <- chEBIToReactomeLowestLevel[grep(ChEBIId, chEBIToReactomeLowestLevel$V1),]
    if (nrow(matchingChEBIToReactome)) {
        TRUE
    } else {
        FALSE
    }
}

#private methods
stringResponseBodyToDataFrame <- function(stringResponseBody, sep = "\t", ...) {
    responseConnection <- textConnection(stringResponseBody)
    dataFrame <- read.table(responseConnection, sep=sep, col.names=c(...),
                                fill=FALSE, strip.white=TRUE)
}
