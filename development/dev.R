###### Depndencies ######
source("http://bioconductor.org/biocLite.R")
# Ubuntu system dependencies
# sudo apt-get install libssl-dev
# sudo apt-get install libxml2-dev

# development dependencies
biocLite(c("devtools", "roxygen2", "knitr", "rmarkdown"))
## Can include important bug fixes.
devtools::install_github("hadley/devtools", force = TRUE)

# build dependencies
biocLite(c("igraph", "pls", "corrplot", "yacca", "FRCC", "RCurl",
           "httr", "RJSONIO", "XML", "S4Vectors", "CCA", "gridExtra",
           "biomaRt", "STRINGdb", "BiocCheck", "S4Vectors",
           "ChemmineOB", "mygene"))

# test dependencies
biocLite("testthat")


###### Vignete Build ######
devtools::use_vignette(name = "OmicsON")
devtools::build_vignettes()
vignette("OmicsON")


###### Unit Tests ######
library(testthat)
devtools::use_testthat()
cat("\014")
devtools::test()


###### Bioconductor - Package Submission ######


###### Integration Tests ######

# CRAN -> R CMD check
library(devtools)
cat("\014")
devtools::check()

# Bioconducto -> R CMD BiocCheck
source("https://bioconductor.org/biocLite.R")
biocLite("BiocCheck")
library(BiocCheck)
cat("\014")
#Do BiocCheck on development version of package, not installed. :P
BiocCheck::BiocCheck(getwd())


###### Writing R documentation files ######
library(roxygen2)
roxygen2::roxygenise()

###### Work Helpers ######
sessionInfo()
R.Version()
cat("\014")
packageVersion("STRINGdb")
find.package("OmicsON")


###### ??? ######
devtools::build(binary = TRUE, args = c('--preclean'))
devtools::build("/home/koralgooll/doktorat/Rpackages/OmicsON/")
devtools::install()
devtools::reload()
install.packages("rstudioapi")
devtools::has_devel()


###### Known Bugs ######

# IMPORTANT note : if some errors restart session or restard IDE!

#Failed to copy the script/BiocCheck script to /usr/lib/R/bin. If you want to be able to
#run 'R CMD BiocCheck' you'll need to copy it yourself to a directory on your PATH, making
#sure it is executable. See the BiocCheck vignette for more information.
# SOLUTION
find.package("BiocCheck")
#Do what is in instruction. Copy that file please. :)


