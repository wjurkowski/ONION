# ======================= rCCA - Regularized Canonical Correlation Analysis  =====================================

# rCCA - rularized extension of the Canonical Correlation Analysis to seek correlations between two data matrices 
# when the number of columns (variables) exceeds the number of rows (observations) 
# When the number of variables is large compared to the number of experimental unitsit is im-
# possible to calculate the inverse of these matrices directly and therefore it is necessary to
# add a multiple of the identity matrix to them. This procedure is known as regularization


# =================================================================================================================
# ------------------------------- FRCCA - fast regularized canonical analysis ------- -----------------------------
# =================================================================================================================


# Fast Regularized Canonical Correlation algorithm described in [Cruz-Cano et al., 2012]. 
# The main idea of the algorithm is using the minimum risk estimators of the correlation 
# matrices described in [Schafer and Strimmer, 2008] during the calculation of the Canonical 
# correlation Structure. It can be considered an extesion of the work for two set of variables 
# (blocks) mentioned in [Tenenhaus and Tenenhaus, 2011] 

options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)

library(stats)
library(corrplot)
library(FRCC)

# --------   IMPORT AND PREPROCESING DATA --------------------------------------
                   
expres <- 
read.table(args[1],
header=TRUE, sep="\t", dec=".", row.names="symbol")
        X <-as.matrix(t(expres))


lipid<-  
read.table(args[2],
header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Y <- as.matrix(lipid)

#--------------- NORMAL DISTRIBUTION TEST (Optionally)
# normal distribution is a prerequisit for application of rCCA
		#X
		lshapX <- apply(X, 2, shapiro.test)
		str(lshapX)
		#Y
		lshapY <- apply(Y, 2, shapiro.test)
		str(lshapY)

#--------------- CORRELATION - reduction number of variables

	#X
		corX <- cor(X)
		corX[upper.tri(corX)] <- 0
		diag(corX) <- 0
		new.X <- X[,!apply(corX,2,function(x) any(x > args[3]))]
		new.X<-new.X[,colSums(is.na(new.X))==0]
		
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
		new.Y<-new.Y[,colSums(is.na(new.Y))==0]
	
			# correlation result for Y
#			ord <- corrMatOrder(corY, order="AOE")
#			M2 <- corY[ord,ord]
#			corrplot.mixed(M2, tl.pos="lt", diag="l")
#			# symnum(corY)

		# X,Y
		corXY <- cor(new.X,new.Y)
		# corrplot(corXY)
		
	tiff("corXY_FRCCA.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
    	corrplot(corXY)
    	dev.off()

#-------------------- DESCRIPTIVE STATISTICS -----------------------------------

		# summaryX <- summary(new.X)
		# summaryY <- summary(new.Y)

		out<-capture.output(summary(new.X))
		cat(out,file="summaryX-FRCCA.txt",sep="\n",append=TRUE)
		out<-capture.output(summary(new.Y))
		cat(out,file="summaryY-FRCCA.txt",sep="\n",append=TRUE)


#------------------- FRCCA CANONICAL ANALYSIS ---------------------------------------- 
	
		
	# my_res <- frcc(X, Y)
	my_res <- frcc(new.X, new.Y)
	# print(my_res) 

	# cor 	- Canonical correlations.
	# p_values 	- The corresponding p-values for the each of the canonical correlations.
	# canonical_weights_X - The canonical weights for the variables of the dataset X.
	# canonical_weights_Y - The canonical weights for the variables of the dataset Y.
	# canonical_factor_loadings_X - The interset canonical factor loadings for the variables of the dataset X.
	# canonical_factor_loadings_Y - The interset canonical factor loadings for the variables of the dataset Y.

# --------	rearrange FRCC for p-value calculation 

	rearrange.frcc(my_res)
		out<-capture.output(my_res )
    		cat(out,file="my_res-FRCCA.txt",sep="\n",append=TRUE)

# This function rearranges the canonical structure according to the canonical correlations from largest to smallest.
# By using the minimum risk estimators of the correlation matrices instead of the sample correlation
# matrices the FRCC algoeithm might disrupt the order of the canonical correlations and hence of
# the canonical structure. This is unacceptable for the algorithm used to calculate the p-values which
# requires the canonical correltions to be ordered in a descending order. This function rearranges the
# canonical structure according to the canonical correlations from largest to smallest.



# --------- graphical presentation
	
	# plot_variables(my_res,1,2, inner_circle_radius = 0.5, text_size = 0.01 )
	
  	tiff("my_res1_FRCCA.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
  	plot_variables(my_res,1,2, inner_circle_radius = 0.5, text_size = 0.01 )
  	dev.off()

	# plot_variables(res.mrcc, i, j, inner_circle_radius = 0.5, text_size = 0.8)
	# Canonical Factor Loadings which will be used as the horizontal axis
	# Canonical Factor Loadings which will be used as the vertical axis
	# Radius of the circle which is used to determine which variables are significant.
	# Only the significant variables will be labled.

	# plot_units(X, Y, my_res, 1, text_size = 0.01, point_size = 1)
	
	tiff("my_res2_FRCCA.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
  	plot_units(new.X, new.Y, my_res, 1, text_size = 0.01, point_size = 1)
  	dev.off()

	# plot_units(X, Y, res.mrcc, i, text_size = 0.8, point_size = 2)
	# i - Canonical Variate which will be used for the axes 
	# (X for horizontal and Y for vertical)
					