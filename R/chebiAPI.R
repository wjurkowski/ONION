library(RCurl)
library(XML)

getChEBIIdsListOfOntologyParents <- function(ChEBIId) {
    headerFields = c(Accept = "text/xml", Accept = "multipart/*", 'Content-Type' = "text/xml; charset=utf-8", SOAPAction = "")
    body = paste0('<soapenv:Envelope
                xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
                xmlns:cheb="http://www.ebi.ac.uk/webservices/chebi">
                <soapenv:Header/>
                <soapenv:Body>
                    <cheb:getOntologyChildren>
                        <cheb:chebiId>', ChEBIId, '</cheb:chebiId>
                    </cheb:getOntologyChildren>
                </soapenv:Body>
            </soapenv:Envelope>', sep = "")
    response <- basicTextGatherer()
    curlPerform(url = "http://www.ebi.ac.uk:80/webservices/chebi/2.0/webservice",
                httpheader = headerFields, postfields = body, writefunction = response$update)
    ontologyParents <- xmlToList(response$value())
    chebiIdsOfontologyParents <- lapply(ontologyParents$Body$getOntologyChildrenResponse$return, function(listElement) {
        strsplit(listElement$chebiId[1], ":")[[1]][2]
    })
}
