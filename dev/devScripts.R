
# Clear R console.
cat("\014")

# origin and upstream integration
# https://help.github.com/articles/configuring-a-remote-for-a-fork/
# https://help.github.com/articles/syncing-a-fork/

# Load all libraries necessary to build package.
setDeveloperEnvironment <- function() {
    install.packages('igraph')
    install.packages('pls')
    install.packages('corrplot')
    install.packages('yacca')
    install.packages('FRCC')
    install.packages('RCurl')
    install.packages('testthat')
    install.packages('httr')
    install.packages('RJSONIO')
    install.packages('XML')
    install.packages('S4Vectors')
    install.packages('CCA')
    install.packages('gridExtra')
    #Since 3.1 bioConductor
    source("http://bioconductor.org/biocLite.R")
    biocLite()
    #install.packages('stats') not available anyweare.
    source('https://bioconductor.org/biocLite.R')
    biocLite('biomaRt')
    biocLite('STRINGdb')
    biocLite('BiocCheck')
    biocLite('S4Vectors')
    biocLite('ChemmineR')
    biocLite('ChemmineOB')
    biocLite("mygene")
}

#TODO: Set up devtools to smart and fast builds.
devtools::install_github("hadley/devtools")
devtools::install_github("Bioconductor-mirror/mygene", force = TRUE)
install.packages("rstudioapi")
library(devtools)
devtools::has_devel()
sessionInfo()

# BioConductor important steps.
# IMPORTANT note : if some errors restart session or restard IDE!
source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
BiocInstaller::useDevel()
BiocInstaller::biocValid()
BiocInstaller::biocLite()
# R CMD BiocCheck
biocLite("BiocCheck")
library(BiocCheck)
#Do BiocCheck on development version of package, not installed. :P
BiocCheck("/home/koralgooll/doktorat/Rpackages/ONION")
#Not this directory:
find.package("ReactomeAPI")
# KNOWN BUG
#Failed to copy the script/BiocCheck script to /usr/lib/R/bin. If you want to be able to
#run 'R CMD BiocCheck' you'll need to copy it yourself to a directory on your PATH, making
#sure it is executable. See the BiocCheck vignette for more information.
# SOLUTION
find.package("BiocCheck")
#Do what is in instruction. Copy that file please. :)

# Build with vignetts.
devtools::use_vignette(name = "ONIONnot")
devtools::build_vignettes()
devtools::build()
devtools::build(binary = TRUE, args = c('--preclean'))
devtools::install()
vignette("ONION")
vignette("ONION")

#TODO: Continous delivery!

#TODO: Set up roxygen2 to autogeneration of manual.

packageVersion("mygene")

devtools::use_testthat()
devtools::test()
