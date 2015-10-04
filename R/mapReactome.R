# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#TODO what arguments? What is the true input? How compicated should be filtering?
mapReactome <- function() {
    (nmLipidomics <- read.table("C:/HOME/ONIONpackage/ONION/R/nm-lipidomics.txt",header=TRUE))
    (chebi <- nmLipidomics[1])
    (chebiFiltered <- chebi["ChEBI"])
    (chebiFiltered <- chebi[chebi$ChEBI > 20000 & chebi$ChEBI < 30000,])
    chebiFiltered <- as.data.frame(chebiFiltered)
    colnames(chebiFiltered) <- c("ChEBI")
    chebiFiltered
}

#TODO Reactome API.
mapReactome()
# http://www.reactome.org/
# http://www.reactome.org/cgi-bin/classbrowser?DB=gk_current&CLASS=Reaction

#install.packages('testthat')
#library(testthat)
