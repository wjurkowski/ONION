#Legend:
#o - operative

#o. Clear R console.
cat("\014")

# origin and upstream integration
# https://help.github.com/articles/configuring-a-remote-for-a-fork/
# https://help.github.com/articles/syncing-a-fork/

#o. Load all libraries necessary to build package.
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
source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script.
biocLite("ChemmineR") # Installs the package.


#TODO: Set up devtools to smart and fast builds.

#TODO: Continous delivery!

#TODO: Set up roxygen2 to autogeneration of manual.
