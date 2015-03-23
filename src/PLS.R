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


# args[1] - path to transcriptional data
# args[2] - path to metabolomics data
# args[3] - number of hidden components  # (standard used ncomp=3)
# args[4] - ścieżka do katalogu, gdzie skrypt powinien wygenerować pliki wynikowe

setwd(commandArgs(TRUE)[4])
options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)

require(pls)
library(pls)

args <- commandArgs(trailingOnly = TRUE)

# ------------------------------- PLS for NUTRIMOUSE --------------------------------------------------

expres <- read.table(args[1], header=TRUE, sep="\t", dec=".", row.names=NULL)
        expres <-as.matrix(t(expres))
        colnames(expres)<-expres[1,]
        expres<-expres[-1,]
        X <- data.matrix(expres)
        class(X) <- "numeric"

lipid <- read.table(args[2], header=TRUE, sep="\t",  dec=".",row.names=NULL)
        lipid <-t(lipid)
        colnames(lipid) <-lipid[1,]
        lipid <- lipid[-1,]
        Y <- as.matrix(lipid)
        class(Y) <- "numeric"



    X <- I(as.matrix(X))
    Y <- I(as.matrix(Y))


    data_train = list(X,Y)


# ---------------- PLS procedure (with model validation) ----------------------------------------

    PLS <- plsr(Y ~ X, data = data_train, ncomp=as.numeric(args[3]), validation="LOO")# "plsr" the same as: "mvr"
#    PLS <- plsr(Y ~ X, data = data_train, ncomp=10, validation="LOO")# "plsr" the same as: "mvr"

    print("past plsr")

    Summary_PLS <- summary(PLS)
        out<-capture.output(Summary_PLS)
        cat(out,file="Summary_PLS.txt",sep="\n",append=TRUE)

        # PL: Sumarycznie przedstawione informacje dotyczące: wymiarów zmiennych w dwóch zbiorach (X i Y), metody integrujacej dane, liczby zmiennych ukrytych i procentu wyjasnionej wariancji.
        # EN: Summary consist of: variable dimension (variables in X and Y sets), fit method, number of components and % variance explained

    print("past summary")

    PLS_t <- loadings(PLS) #  współczynniki dla macierzy utajonych względem X   # t
        out<-capture.output(PLS_t)
        cat(out,file="PLS_t.txt",sep="\n",append=TRUE)

        # PL: Ładunki czynnikowe (korelacje między zmiennymi utajonymi a zmiennymi w każdym zbiorze)
        # EN: Loadings (the correlations between the latent variables and each variables in the dataset )

    print("past loadings")

    PLS_loadings <- loading.weights(PLS)
        out<-capture.output(PLS_loadings)
        cat(out,file="PLS_loadings.txt",sep="\n",append=TRUE)

        # PL: Macierz wag czynnikowych (dla zbioru X), gdzie każdy wiersz zawiera współczynniki definiujące liniowa kombinację zmiennych ukrytych (LV)(bezwzględne wartości wag odzwierciedlają udział zmiennych w modelu)
        # EN: matrix of predictor loadings (for X set), where each row contains coefficients that define a linear combination of PLS components (LV) that approximate the original predictor variables.


    print("past weights")

    PLS_T <- scores(PLS) # macierz T
        out<-capture.output(PLS_T)
        cat(out,file="PLS_T.txt",sep="\n",append=TRUE)

        # PL: Macierz wartości czynnikowych dla odpowiedniej macierzy wag zbioru X (T=XW) (macierz ortogonalna z wierszami odpowiadającymi obserwacjom i zmiennymi ukrytymi w kolumnach)
        # EN: The scores matrix reflect weight matrix W for X such that T=XW (orthonormal matrix with rows corresponding to observations, columns to components)

    print("past scores")

    PLS_U <- Yscores(PLS) # macierz U
        out<-capture.output(PLS_U)
        cat(out,file="PLS_U.txt",sep="\n",append=TRUE)

        # PL: Macierz wartości czynnikowych dla odpowiedniej macierzy wag zbioru Y (macierz ortogonalna z wierszami odpowiadającymi obserwacjom i zmiennymi ukrytymi w kolumnach)
        # EN: The scores matrix reflect weight matrix for Y set of data (orthonormal matrix with rows corresponding to observations, columns to components)

    print("past yscores")

    PLS_Y <- Yloadings(PLS) # współczynniki dla macierzy utajonych względem Y   # u
        out<-capture.output(PLS_Y)
        cat(out,file="PLS_Y.txt",sep="\n",append=TRUE)

        # PL: macierz wagi czynnikowych zbioru Y, gdzie każdy wiersz zawiera współczynniki definiujące liniowa kombinację zmiennych ukrytych (LV)(bezwzględne wartości wag odzwierciedlają udział zmiennych w modelu)
        # EN: matrix of response loadings (for Y set), where each row of YLOADINGS contains coefficients that define a linear combination of PLS components that approximate the original response variables

    print("past yloadings")

#    PLS_B <- coef(PLS, ncom = 1:3, legendpos = "bottomright",data = data_train)   # B
    PLS_B <- coef(PLS, ncom = 1:3, data = data_train)   # B
        out<-capture.output(PLS_B)
        cat(out,file="PLS_B.txt",sep="\n",append=TRUE)

        # PL: Macierz współczynników regresji (B)
        # EN: Regression coefficient matrix (B)

    print("past coef")


# ----------- graphical presentation -------------------------------------------------------------

# plot(x, plottype = c("prediction", "validation", "coefficients",
# "scores", "loadings", "biplot", "correlation"), ...)


#------TIFF --------------------

        # tiff("PLS_RMSEP.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
            # plot(PLS, "validation", val.type = "RMSEP") # error should be low -> good for prediction
        # dev.off()

            # validationplot(PLS, estimate = "all")

        # tiff("PLS_loadings.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
            # par(mfrow = c(2,2))
            # biplot(PLS, which = "x") # Default
            # biplot(PLS, which = "y")
            # biplot(PLS, which = "scores")
            # biplot(PLS, which = "loadings")
        #dev.off()

        # tiff("PLS_loadings_variables.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
            # loadingplot(PLS, comps = 1:3, legendpos = "topright") # With legend
        # dev.off()

        # tiff("PLS_loadings_scaterrplot.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
            # loadingplot(PLS, comps = 1:3, scatter = TRUE) # Plot as scatterplots
        # dev.off()

        # tiff("PLS_coefplot.tiff", res=600, compression = "lzw", height=5, width=5, units="in")
            # coefplot(PLS, ncom = 1:3, legendpos = "bottomright",data = data_train)
        # dev.off()


#------PNG --------------------

#    png("PLS_RMSEP.png", width = 640, height = 480)
#        plot(PLS, "validation", val.type = "RMSEP") # error should be low -> good for prediction
#    dev.off()

        # PL Wykresy RMSEP (an. root mean squared error of prediction). Im mniejsza wartość RMSEP tym lepsze dopasowanie modelu.
        # EN Plots of root mean squared error of prediction (RMSEP). The lowest value of RMSEP the better fit.

    png("PLS_loadings.png", width = 640, height = 480)
        par(mfrow = c(2,2))
        biplot(PLS, which = "x") # Default
        biplot(PLS, which = "y")
        biplot(PLS, which = "scores")
        biplot(PLS, which = "loadings")
    dev.off()

        # PL Wykresy: współczynników X i ładunków czynnikowych zbioru X, współczynników Y i ładunków czynnikowych zbioru Y, współczynników X i współczynników Y, ładunków czynnikowych zbioru X i współczynników X i ładunków czynnikowych zbioru Y,
#            /Współczynniki X: Macierz wartości czynnikowych dla odpowiedniej macierzy wag zbioru X (T=XW) (macierz ortogonalna z wierszami odpowiadającymi obserwacjom i zmiennymi ukrytymi w kolumnach)
#            /Ładunki czynnikowe: korelacje między zmiennymi utajonymi a zmiennymi w każdym zbiorze/
#            /współczynniki Y: Macierz wartości czynnikowych dla odpowiedniej macierzy wag zbioru Y (macierz ortogonalna z wierszami odpowiadającymi obserwacjom i zmiennymi ukrytymi w kolumnach)/
#            /ładunki czynnikowe zbioru Y: macierz wag czynnikowych zbioru Y, gdzie każdy wiersz zawiera współczynniki definiujące liniowa kombinację zmiennych ukrytych (LV)(bezwzględne wartości wag odzwierciedlają udział zmiennych w modelu)/

        # EN Plots of: X scores and X loadings, Y scores and Y loadings, X scores and Y scores, X loadings and Y loadings
#            /X scores: The scores matrix reflect weight matrix W for X such that T=XW (orthonormal matrix with rows corresponding to observations, columns to components)/
#            /Loadings: the correlations between the latent variables and each variables in the dataset/
#            /Y scores: The scores matrix reflect weight matrix for Y set of data (orthonormal matrix with rows corresponding to observations, columns to components)/
#            /Y loadings: matrix of response loadings (for Y set), where each row of YLOADINGS contains coefficients that define a linear combination of PLS components that approximate the original response variables/


    png("PLS_loadings_variables.png", width = 640, height = 480)
    loadingplot(PLS, comps = 1:3, legendpos = "topright") # With legend
    dev.off()

        # PL: Wykres Ładunków czynnikowych dla zbioru X z informacją o procencie wyjasnonej wariancji przez poszczególne zmienne ukryte
        # EN: Plot of Loadings: the correlations between the latent variables and each variables in the X dataset. Percent of variances explained for three latent variables is shown.

    png("PLS_loadings_scaterrplot.png", width = 640, height = 480)
    loadingplot(PLS, comps = 1:3, scatter = TRUE) # Plot as scatterplots
    dev.off()

        # PL: Wykres ładunków czynnikowych dla kazdej zmiennej ukrytej
        # EN: Scatterplot of loadings for each latent variables

    png("PLS_coefplot.png", width = 640, height = 480)
    coefplot(PLS, ncom = 1:3, legendpos = "bottomright",data = data_train)
    dev.off()

        # PL: Wykresy współczynników regresji względem każdej zmiennej ze zbioru Y
        # EN: Plots of correlation coefficient for each of variables in Y set.

    # coefplot(PLS, ncomp = 1:3, separate = TRUE)
    # biplot(PLS) # equivlent: plot(PLS, plottype = "biplot")
    # str(PLS)
