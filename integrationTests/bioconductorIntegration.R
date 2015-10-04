#o. Clear R console.
cat("\014")

library(BiocCheck)
pathToInstallableONIONPackage <- file.path("C:/HOME/doktorat","ONION")
BiocCheck(pathToInstallableONIONPackage)

pathToONIONLibrary <- file.path("C:/HOME/R-3.2.2/library","ONION")