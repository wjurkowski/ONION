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
