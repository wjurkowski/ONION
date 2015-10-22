source("R/reactomeAPI.R")
source("R/chebiAPI.R")

#Original pseudo code API (OPCAPI):
#
#data.frame clusterSmallMolecules( file ) {
#    //file - nm-lipidomics.txt (tylko male czasteczki)
#    //K klastrów (Rodzice Chebi)
#    //  H1    H2
#    //1 123   ABC
#    //2 231   DEF
#    //NIE do pliku, jako obiekt.
#    //output - fa_mapping.txt (dzieci, rodzice) w ID z Chebi i chebi parent by ontology
#}

clusterSmallMolecules <- function(pathToFile, header=TRUE) {
    smallMolecules <- read.table(pathToFile,header)
    smallMoleculesParents <- lapply(as.list(smallMolecules$ChEBI), function(listElement) {
        if (checkIfPathwayIdExistsForChEBIId(listElement)) {
            listElement
        } else {
            listOfChEBIIdsOfOntologyParents <- getListOfChEBIIdsOfOntologyParents(listElement)
            bestParentId <- NA
            for (i in 1:length(listOfChEBIIdsOfOntologyParents)) {
                j <- listOfChEBIIdsOfOntologyParents[[i]]
                if (checkIfPathwayIdExistsForChEBIId(j)) {
                    bestParentId <- j
                    break
                }
            }
            bestParentId
        }
    })
    chEBIIdsToChEBIParentsIDs <- data.frame(ChEBI=as.character(as.list(smallMolecules$ChEBI)), bestChEBIParents=as.character(smallMoleculesParents))
    chEBIIdsToChEBIParentsIDs
}

#Go to chebi, and find at least good identifier in ontology.
#Go to http://www.reactome.org/download/current/ChEBI2Reactome.txt and get Pathways Ids related to CHebi.
#http://www.reactome.org/cgi-bin/instancebrowser?DB=gk_current&ID=6816043&  -  here you can obtain Genes, Chebis and etc.


#OPCAPI:
#
#list[list] mapReactomePathways <- function ( data.frame fa_mapping.txt, ) {
#    //input - fa_mapping.txt
#    //all mapReactome.R file
#    //MY: Find not duplicated PAthways
#    //TODO: find Gens.
#    //smallMolecules[Genes], NOW:smallMolecules[Pathways]
#    //output genesSymbolList smallMolecules[Genes] K[x1...xn]
#}

mapReactomePathways <- function(chEBIToParentsIdDataFrame) {
    parentIdsList <- as.list(as.vector(chEBIToParentsIdDataFrame$bestChEBIParents))
    parentIdsToPathways <- lapply(parentIdsList, function(listElement) {
        if (listElement != "NA") {
            listOfUsablePathwaysIds <- getUsablePathwaysIdsForChEBI(listElement)
        } else {
            NA
        }
    })
    parentIdsToPathways
#    TODO: Map PAthways to Genes. But how?
#    15378	R-ATH-211976	http://www.reactome.org/PathwayBrowser/#R-ATH-211976
#    Endogenous sterols	IEA	Arabidopsis thaliana
#    AT1G17060.1	R-ATH-211976	http://www.reactome.org/PathwayBrowser/#R-ATH-211976
#    Endogenous sterols	IEA	Arabidopsis thaliana

#    candidateSet 2975806 //Reactome - Set
#    entityWithAccessionedSequence 193137 //ENSEMBL, UniProt - Protein
#    simpleEntity 113533 //ChEBI - Chemical Compound
#    complex 5580259 //Reactome - Complex
}
