#o. Clear R console.
cat("\014")

#Run R CMD CHECK.
install.packages('devtools')
library('devtools')
check(document = FALSE)

#Create and Run tests.
library(devtools)
#It will ask you about test path creation to testthat.R
#test()

#Set proper work directory.



#Analyze to understend lists.
dane <- getChEBIIdsListOfOntologyParents(28842)
dane

library(ONION)
ONION::mapReactome()
cat("\014")


#TODO bug with ChEBI names.
l <- list(aaa=c("as as as"), aaa=c("aa bb cc"))
l$aaa
