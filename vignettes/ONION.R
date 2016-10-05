## ----global_options, include=TRUE----------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)

## ---- results = 'asis'---------------------------------------------------
pathToFileWithChebiIds <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")
M <- read.table(pathToFileWithChebiIds, header = TRUE)
pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
T <- read.table(pathToExampleFileWithXData, header = TRUE)
tab <- head(M,5)
knitr::kable(tab[1:7], caption = "Metabolomics data")

## ---- echo=TRUE, eval=TRUE-----------------------------------------------
library(ONION)

## ---- echo=TRUE, results='hide'------------------------------------------
find.package("ONION")

## ---- echo=TRUE, results='asis'------------------------------------------
clusteredSmallMolecules <- ONION::clusterUsingOntology(chebiIdsDataFrame = M,
                                                           rootColumnName = "ChEBI",
                                                           ontologyRepresentatnion = ONION::firstExistsInReactomeChebiOntology)
    
knitr::kable(head(clusteredSmallMolecules, 6))

## ---- echo=FALSE, results='asis'-----------------------------------------
mergedSmallMolecules <- ONION::mergeChEBIOntologyWithChildFavoring(clusteredSmallMolecules, rootColumnName = 'root')
knitr::kable(head(mergedSmallMolecules, 6))
#save mapped IDs
merged <- mergedSmallMolecules$ontologyId[mergedSmallMolecules$ontologyId > 0]
write.table(merged, "metabolites_in_Reactome.txt", quote = F, row.names = F, col.names = F)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(data.frame(mergedSmallMolecules[1:4, c("ontologyId"), drop = FALSE], XYZ = c("1","2","3","4")))

## ---- echo=TRUE, results='hide'------------------------------------------
chebiIdsToReactomePathways <- invisible(ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[, c("ontologyId"), drop = FALSE], organismTaxonomyId = '9606', idsColumnName = "ontologyId", rootColumnName = NULL))
chebiIdsToReactomePathwaysWithRoot <- invisible(ONION::mapReactomePathwaysUnderOrganism(chebiOntologyIds = mergedSmallMolecules[, c("ontologyId", "root"), drop = FALSE], organismTaxonomyId = '9606', idsColumnName = "ontologyId", rootColumnName = "root"))

## ---- echo=FALSE, results='asis'-----------------------------------------
head(chebiIdsToReactomePathwaysWithRoot$ontologyId, 1)
# TODO : Convert root to character in code of ONION.
head(as.character(chebiIdsToReactomePathwaysWithRoot$root), 1)
chebiIdsToReactomePathwaysWithRoot$ensembleIds[[1]]
chebiIdsToReactomePathwaysWithRoot$uniProtIds[[1]]
chebiIdsToReactomePathwaysWithRoot$reactomeIds[[1]]
chebiIdsToReactomePathwaysWithRoot$genesSymbolsFromEnsemble[[1]]
chebiIdsToReactomePathwaysWithRoot$genesSymbolsFromUniProt[[1]]
# knitr::kable(head(chebiIdsToReactomePathways, 1))

## ---- echo=FALSE, results='hide'-----------------------------------------
chebiIdsToReactomePathwaysAndToStringNeighbours <- invisible(ONION::getStringNeighbours(chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",],
                                                                                  stringOrganismId = 9606,
                                                                                  stringDbVersion = "10",
                                                                                  idsColumnName = 'ontologyId',
                                                                                  rootColumnName = NULL,
                                                                                  listOfEnsembleIdColumnName = 'ensembleIds'))

## ---- echo=FALSE, results='asis'-----------------------------------------
# knitr::kable(chebiIdsToReactomePathwaysAndToStringNeighbours[chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId == "CHEBI:15756",])
chebiIdsToReactomePathwaysAndToStringNeighbours$ontologyId
chebiIdsToReactomePathwaysAndToStringNeighbours$ensembleIds[[1]]
chebiIdsToReactomePathwaysAndToStringNeighbours$stringIds[[1]]
chebiIdsToReactomePathwaysAndToStringNeighbours$stringGenesSymbols[[1]]

## ---- echo=FALSE, results='asis'-----------------------------------------
gmtGroupsFilePath <- paste(find.package("ONION"),"/example/nm-groups.txt", sep = "")
groups <- ONION::groupUsingUserDefinition(pathToFileWithGroupDefinition = gmtGroupsFilePath)
knitr::kable(groups)

## ---- echo=FALSE, results='hide', eval=FALSE-----------------------------
#  #select small molecules
#  lip1 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:27432",]$root
#  lip2 <- mergedSmallMolecules[mergedSmallMolecules$root == "CHEBI:73705",]$root
#  joinLip <- c(as.character(lip1), as.character(lip2))
#  
#  #use Reactome genes mapped to selected small molecules
#  reactomeTrans1 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:15756",]$genesSymbolsFromEnsemble[[1]]
#  reactomeTrans2 <- chebiIdsToReactomePathways[chebiIdsToReactomePathways$ontologyId == "CHEBI:16015",]$genesSymbolsFromEnsemble[[1]]
#  joinRecatomeTrans <- c(reactomeTrans1, reactomeTrans2)[!duplicated(c(reactomeTrans1, reactomeTrans2))]

## ---- echo=FALSE, results='hide', eval=FALSE-----------------------------
#  functionalInteractions <- ONION::createFunctionalInteractionsDataFrame(chebiIdsToReactomePathways)
#  knitr::kable(head(functionalInteractions, 6))

## ---- echo=FALSE, results='hide', eval=FALSE-----------------------------
#      pathToExampleFileWithXData <- paste(find.package("ONION"),"/example/nm-transcriptomics.txt", sep = "")
#      pathToExampleFileWithYData <- paste(find.package("ONION"),"/example/nm-lipidomics.txt", sep = "")
#  
#      XDF <- read.table(pathToExampleFileWithXData, header = TRUE);
#      YDF <- read.table(pathToExampleFileWithYData, header = TRUE);
#  
#      ccaResults1 <- ONION::makeCanonicalCorrelationAnalysis(
#          xNamesVector = joinRecatomeTrans,
#          yNamesVector = joinLip,
#              XDataFrame = XDF,
#              YDataFrame = YDF)

