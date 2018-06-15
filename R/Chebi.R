
getChEBIOntologyChildren <- function(ChEBIId) {
    headerFields = c(Accept = "text/xml", Accept = "multipart/*", 'Content-Type' = "text/xml; charset=utf-8", SOAPAction = "")
    body = paste0('<soapenv:Envelope
                  xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
                  xmlns:cheb="https://www.ebi.ac.uk/webservices/chebi">
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
    ontologyChildren <- xmlToList(response$value())
    ChEBIOntologhyChildren <- ldply(ontologyChildren$Body$getOntologyChildrenResponse$return, function(listElement){
        listElement$chebiId <- strsplit(as.character(listElement$chebiId), ":")[[1]][2]
        df <- data.frame(listElement)
        df
    })
    ChEBIOntologhyChildren
}

getChEBIOntologyParents <- function(ChEBIId) {
    headerFields = c(Accept = "text/xml", Accept = "multipart/*", 'Content-Type' = "text/xml; charset=utf-8", SOAPAction = "")
    body = paste0('<soapenv:Envelope
                  xmlns:soapenv="http://schemas.xmlsoap.org/soap/envelope/"
                  xmlns:cheb="https://www.ebi.ac.uk/webservices/chebi">
                  <soapenv:Header/>
                  <soapenv:Body>
                    <cheb:getOntologyParents>
                        <cheb:chebiId>', ChEBIId, '</cheb:chebiId>
                    </cheb:getOntologyParents>
                  </soapenv:Body>
                  </soapenv:Envelope>', sep = "")
    response <- basicTextGatherer()
    curlPerform(url = "http://www.ebi.ac.uk:80/webservices/chebi/2.0/webservice",
                httpheader = headerFields, postfields = body, writefunction = response$update)
    ontologyParents <- xmlToList(response$value())
    ChEBIOntologhyParents <- ldply(ontologyParents$Body$getOntologyParentsResponse$return, function(listElement){
        listElement$chebiId <- strsplit(as.character(listElement$chebiId), ":")[[1]][2]
        df <- data.frame(listElement)
        df
    })
    ChEBIOntologhyParents
}
