# =================================================================================================================
# ======================= CCA - canonical analysis  ===============================================================
# =================================================================================================================

# Canonical analysis is a multivariate technique which is concerned with determining the relationships 
# between groups of variables in a data set. The data set is split into two groups, let's call these 
# groups X and Y, based on some common characteristics. The purpose of Canonical analysis is then 
# to find the relationship between X and Y, i.e. can some form of X represent Y. It works by finding 
# the linear combination of X variables, i.e. X1, X2 etc., and linear combination of Y variables, 
# i.e. Y1, Y2 etc., which are most highly correlated. This combination is known as the "first canonical 
# variates" which are usually denoted U1 and V1, with the pair of U1 and V1 being called a "canonical function". 
# The next canonical functions, U2 and V2 are then restricted so that they are uncorrelated with U1 and V1.
# Everything is scaled so that the variance equals 1. 

# more about yacca option http://cran.r-project.org/web/packages/yacca/yacca.pdf

options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)

library(stats)
library(corrplot)
library(yacca)


# ------------------------------- Canonical analysis tutorial example -----------------------------------------------------------

####  Example parallels the R builtin cancor example
	# data(LifeCycleSavings)
	# pop <- LifeCycleSavings[, 2:3]
	# oec <- LifeCycleSavings[, -(2:3)]
	# cca.fit <- cca(pop, oec)

#### View the results
	# cca.fit
	# summary(cca.fit)
	# plot(cca.fit)


# --------   IMPORT DATA  -----------------------------------------------

expres <- read.table(args[1],header=TRUE, sep="\t", dec=".", row.names="symbol")
        X <-as.matrix(t(expres))


lipid<-  read.table(args[2],header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Y <- as.matrix(lipid)


#!!!krok wymaga automatyzacji w celu wybrania zmiennych ktorych rozklad jest zgodny z normalnym
	# normally distribution test
		#X
		lshapX <- apply(X, 2, shapiro.test)
		str(lshapX)
		#Y
		lshapY <- apply(Y, 2, shapiro.test)
		str(lshapY)


	# correlation. Reduction number of variables on the basis of correlation results
	#X
		corX <- cor(X)
		corX[upper.tri(corX)] <- 0
		diag(corX) <- 0
		new.X <- X[,!apply(corX,2,function(x) any(x > args[3]))]
		# head(new.X)
		# dim(new.X)
		# dim(X)
		
			# correlation result for X
#			ord <- corrMatOrder(corX, order="AOE")
#			M2 <- corX[ord,ord]
#			corrplot.mixed(M2, tl.pos="lt", diag="l")
#			symnum(corX)

	#Y
		corY <- cor(Y)
		corY[upper.tri(corY)] <- 0
		diag(corY) <- 0
		new.Y <- Y[,!apply(corY,2,function(x) any(x > args[4]))]
		# head(new.Y)
		# str(new.Y)
		# str(Y)

			# correlation result for Y
#			ord <- corrMatOrder(corY, order="AOE")
#			M2 <- corY[ord,ord]
#			corrplot.mixed(M2, tl.pos="lt", diag="l")
#			# symnum(corY)

		# X,Y
		corXY <- cor(new.X,new.Y)
		corrplot(corXY)
		
		png("corXY.png", width = 640, height = 480) 
		corrplot(corXY)
		dev.off()

	# descriptive statistics
		# summaryX <- summary(new.X)
		# summaryY <- summary(new.Y)

		out<-capture.output(summary(new.X))
		cat(out,file="summaryX-YACCA.txt",sep="\n",append=TRUE)
		out<-capture.output(summary(new.Y))
		cat(out,file="summaryY-YACCA.txt",sep="\n",append=TRUE)


# ------------------------------
		
	cca.fit <- cca(new.X, new.Y)
	out<-capture.output(cca.fit)
	cat(out,file="cca.fit-YACCA.txt",sep="\n",append=TRUE)


#	cca.fit <- cca(X, Y)  
	
# results
	cca.fit 
	# => canonical correlations (CV);
	# => X (or Y) Coefficient = Coefficients for the x (or y) variables 
	# on each canonical variate (canonical weight) (PL: wagi okreœlaj¹ce 
	# zmienne kanonicze. Im wiêksza bezwzglêdna wartoœæ wagi wym wiekszy 
	# wk³ad danej zmiennej do modelu)
	# => structural correlation loadings (PL: kanoniczne ³adunki czynnikowe) 
	# Aggregate Redundancy Coefficients (coefficient accounted for by X (or Y) 
	# variables, through all canonical variates (PL: ile przecietnej wariancji w jenym zbiorze jest wyjaœnione przez
	# dana zmienna kanonicza w drugim zbiorze) 

	summary(cca.fit) # wymaga wybrania interesuj¹cych danych

	png("corXY.png", width = 640, height = 480) 
	plot(cca.fit)
	dev.off()

# test for significance canonical correlations
	# Several related tests have been proposed for the evaluation of canonical correlations (including
	# Bartlett’s Chi-squared test, which is computed by default within cca). This function employs Rao’s
	# statistic (related to Wilks’ Lambda) as the basis for an F test of each coefficient (and all others in
	# ascending sequence) against the hypothesis that the associated population correlations are zero
	
	# F.test.cca(cca.fit) 
	out<-capture.output(F.test.cca(cca.fit))
	cat(out,file="F_test_cca.txt",sep="\n",append=TRUE)

# helio plot - show loadings on first canonical variate
	# Helio plots display data in radial bars, with larger values pointing outward from a base reference
	# circle and smaller (more negative) values pointing inward). Such plots are well-suited to the display
	# of multivariate information with several groups of variables, as with canonical correlation analysis
	
	# helio.plot(cca.fit, x.name="GENES",y.name="LIPIDS")
	png("helio_plot_YACCA.png", width = 640, height = 480) 
	helio.plot(cca.fit, x.name="GENES",y.name="LIPIDS")
	dev.off()

# helio plot - show variances on second canonical variate
	helio.plot(cca.fit, cv=2, x.name="GENES",y.name="LIPIDS", type="variance")


	
