<<<<<<< HEAD
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

# ------------------------------- IMPORT DATA --------------------------------------------------# -------------------- Import data

expres <- 
read.table("args[1]",
header=TRUE, sep="\t", dec=".", row.names="symbol")
        X <-I(as.matrix(t(expres)))
	

lipid<-  
read.table("args[2]",
header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Y <- I(as.matrix(lipid))
	
	
	# expres <- read.csv("/project/integromika/data/nutrimouse/grupa_nutrimouse/grupaexp_nutrimouse.csv",header=TRUE, sep="\t", dec=".")	
	#	expres <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseGENE.csv", header=TRUE, sep=";", dec=".", row.names = NULL)
	#	lipid  <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseLIPID.csv", header=TRUE, sep=";", dec=".",row.names = NULL)			


	data_train = list(X,Y)

# ---------------- PLS procedure ------------------------------------------------------------------

	PLS<- plsr(Y ~ X, data = data_train, validation = "CV", ncomp=3)
	summary(PLS)

# ----------- model validation

# plot(x, plottype = c("prediction", "validation", "coefficients",
# "scores", "loadings", "biplot", "correlation"), ...)

plot(PLS, "val", val.type = "MSEP", estimate = "CV")
validationplot(PLS, estimate = "all")
plot(PLS, "validation", val.type = "RMSEP") # error should be low -> good for prediction


# ----------- canonical values

	loadings(PLS) #  wspó³czynniki dla macierzy utajonych wzglêdem X   # t
	loading.weights(PLS)
	scores(PLS) # macierz T
	Yscores(PLS) # macierz U
	Yloadings(PLS) # wspó³czynniki dla macierzy utajonych wzglêdem Y   # u
	PLS$projection # P ?????

# ------------ graphical presentation:
	biplot(PLS) # equivlent: plot(PLS, plottype = "biplot")
	# str(PLS)
	
	wsp.cor <- coef(PLS, ncom = 1:3, legendpos = "bottomright",data = data_test)   # B
	coefplot(PLS, ncom = 1:3, legendpos = "bottomright",data = data_test)
	coefplot(PLS, ncomp = 1:3, separate = TRUE)

	loadingplot(PLS, comps = 1:3, legendpos = "topright") # With legend
	loadingplot(PLS, comps = 1:3, scatter = TRUE) # Plot as scatterplots
	

# The four combinations of x and y points:
par(mfrow = c(2,2))
	biplot(PLS, which = "x") # Default
	biplot(PLS, which = "y")
	biplot(PLS, which = "scores")
	biplot(PLS, which = "loadings")





=======
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

# ------------------------------- IMPORT DATA --------------------------------------------------# -------------------- Import data

expres <- 
read.table("args[1]",
header=TRUE, sep="\t", dec=".", row.names="symbol")
        X <-I(as.matrix(t(expres)))
	

lipid<-  
read.table("args[2]",
header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Y <- I(as.matrix(lipid))
	
	
	# expres <- read.csv("/project/integromika/data/nutrimouse/grupa_nutrimouse/grupaexp_nutrimouse.csv",header=TRUE, sep="\t", dec=".")	
	#	expres <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseGENE.csv", header=TRUE, sep=";", dec=".", row.names = NULL)
	#	lipid  <-  read.csv("/project/integromika/data/nutrimouse/nutrimouseLIPID.csv", header=TRUE, sep=";", dec=".",row.names = NULL)			


	data_train = list(X,Y)

# ---------------- PLS procedure ------------------------------------------------------------------

	PLS<- plsr(Y ~ X, data = data_train, validation = "CV", ncomp=3)
	summary(PLS)

# ----------- model validation

# plot(x, plottype = c("prediction", "validation", "coefficients",
# "scores", "loadings", "biplot", "correlation"), ...)

plot(PLS, "val", val.type = "MSEP", estimate = "CV")
validationplot(PLS, estimate = "all")
plot(PLS, "validation", val.type = "RMSEP") # error should be low -> good for prediction


# ----------- canonical values

	loadings(PLS) #  wspó³czynniki dla macierzy utajonych wzglêdem X   # t
	loading.weights(PLS)
	scores(PLS) # macierz T
	Yscores(PLS) # macierz U
	Yloadings(PLS) # wspó³czynniki dla macierzy utajonych wzglêdem Y   # u
	PLS$projection # P ?????

# ------------ graphical presentation:
	biplot(PLS) # equivlent: plot(PLS, plottype = "biplot")
	# str(PLS)
	
	wsp.cor <- coef(PLS, ncom = 1:3, legendpos = "bottomright",data = data_test)   # B
	coefplot(PLS, ncom = 1:3, legendpos = "bottomright",data = data_test)
	coefplot(PLS, ncomp = 1:3, separate = TRUE)

	loadingplot(PLS, comps = 1:3, legendpos = "topright") # With legend
	loadingplot(PLS, comps = 1:3, scatter = TRUE) # Plot as scatterplots
	

# The four combinations of x and y points:
par(mfrow = c(2,2))
	biplot(PLS, which = "x") # Default
	biplot(PLS, which = "y")
	biplot(PLS, which = "scores")
	biplot(PLS, which = "loadings")





>>>>>>> b2f19d081b947847f2fe095549537237be7fd748
