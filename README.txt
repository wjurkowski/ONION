Download OmicsON_0.0.1.tar.gz from builds directory please.

Before installation, remember that OmicsON depends from two packages not from CRAN:
ChebiAPI, ReactomeAPI

Then use command below to install package (remeber to change file path):
install.packages("/home/koralgooll/doktorat/Rpackages/OmicsON_0.0.1.tar.gz", repos = NULL, type="source")

I case of some problem, at first install packages from DESCRIPTION file please:
igraph, ChebiAPI, plyr, ReactomeAPI
RCurl, RJSONIO, XML, testthat, pls, yacca, FRCC, httr, biomaRt, STRINGdb, BiocCheck, CCA, gridExtra, mygene
To installation you can use setDeveloperEnvironment() function from dev/devScripts.R, but remember that you need to look at console and type a and y on quesrions.
