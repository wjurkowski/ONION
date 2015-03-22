# =================================================================================================================
# ======================= PLS - Partial Least Squares Regression  =================================================
# =================================================================================================================

# Partial least squares regression (PLS regression) is a statistical method that bears some relation 
# to principal components regression; instead of finding hyperplanes of minimum variance between 
# the response and independent variables, it finds a linear regression model by projecting the predicted 
# variables and the observable variables to a new space. 
# PLS regression is today most widely used in chemometrics and related areas. 
# It is also used in bioinformatics, sensometrics, neuroscience and anthropology.
# http://cran.r-project.org/web/packages/pls/index.html

ptions(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)

	library(pls)

# ------------------------------- PLS for NUTRIMOUSE --------------------------------------------------# -------------------- Import data

expres <- 
read.table("args[1]",
header=TRUE, sep="\t", dec=".", row.names="symbol")
        Xt <-I(as.matrix(t(expres)))


lipid<-  
read.table("args[2]",
header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Yt <- I(as.matrix(lipid))
	
	
expres <- 
read.table("args[1]",
header=TRUE, sep="\t", dec=".", row.names="symbol")
        Xp <-as.matrix(t(expres))


lipid<-  
read.table("args[2]",
header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Yp <- as.matrix(lipid)
	
#	expres <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseGENE.csv", header=TRUE, sep=";", dec=".", row.names = NULL)
#	lipid  <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseLIPID.csv", header=TRUE, sep=";", dec=".",row.names = NULL)			
#	X <- I(as.matrix(expres[,-1]))
#	Y<- I(as.matrix(lipid[,-1]))

	data_train = list(Xt,Yt)
	data_test =  list(Xp,Yp)


# --------------------------- PLS for prediction -------------------------------------------------------

	PLS.mvr <- mvr(Yt ~ Xt, data = data_train, validation = "CV", ncomp=3)
	
	summary(PLS.mvr)

## Predicted responses for models with 1, 2, 3 and 4 components
	pred.resp <- predict(PLS.mvr, ncomp = 1:4, newdata = data_test)

## Predicted scores for components
	predict(PLS.mvr, comps = 1:3, type = "scores", newdata = data_test)
	
## Both cross-validated and test set predictions:
	predplot(PLS.mvr , ncomp = 1:3, which = c("validation", "test"),
	newdata = data_train)
	
	

