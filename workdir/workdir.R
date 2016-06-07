
cat('\014')



xNamesVector = c("15756")
yNamesVector = c("IDI1", "ACACB", "FDFT1", "SQLE", "SC5D", "MVK", "HMGCS1", "HMGCR", "GPAM", "TM7SF2", "GGPS1", "LSS", "FDPS", "PMVK", "MVD", "FASN", "ELOVL6", "DHCR7", "ABCB4", "ALAS1", "ME1", "TRIB3", "SULT2A1", "CPT1A", "APOA5", "FHL2", "ACADM", "APOA1", "CTGF", "G0S2", "SLC27A1", "HMGCS2", "GRHL1", "CD36", "AGT", "GLIPR1", "CYP1A1", "RGL1", "TNFRSF21", "TIAM2", "PLIN2", "ANKRD1", "FADS1", "ACSL1", "CPT2", "APOA2", "ACOX1", "FABP1", "ABCA1", "PEX11A", "ANGPTL4", "CYP7A1", "NPAS2", "PPARA", "CYP4A11", "TXNRD1", "UGT1A9", "KLF5", "PCK1", "PPARG", "EBF1", "PLIN1", "FABP4", "CEBPB", "LEP", "LPL", "ADIPOQ", "SLC2A4", "CEBPD", "CEBPA", "GCG", "TFB1M", "POLRMT", "SSBP1", "NRF1", "C10orf2", "TFAM", "USP46", "ATP5B", "MTERF1", "SIRT3", "GABPA", "TFB2M", "CYCS", "ESRRA", "POLG2", "CRY1", "F7", "RORA", "AVP", "DBP", "NAMPT", "SERPINE1", "CRY2", "BHLHE41", "NR1D1", "PER2", "BHLHE40", "NOCT", "PER1", "GIP", "SREBF1", "ARNTL", "CLOCK", "NPPA")
yNamesVector = chebiIdsToReactomePathways[chebiIdsToReactomePathways$chebiId == 27432,]$gensSymbols[[1]]

pathToFileWithXData = "/home/koralgooll/doktorat/Rpackages/ONION/example/nm-transcriptomics.txt"
pathToFileWithYData = "/home/koralgooll/doktorat/Rpackages/ONION/example/nm-lipidomics-valid.txt"

# Where XData = transcriptomicsData and YData = lipidomicsData.
XData <- readWithoutDuplicates(pathToFileWithXData)
YData <- readWithoutDuplicates(pathToFileWithYData)
str(XData)

transposedXData <- as.data.frame(t(XData))
str(transposedXData)
transposedYData <- as.data.frame(t(YData))

interX <- intersect(colnames(transposedXData), xNamesVector)
interY <- intersect(colnames(transposedYData), joinLip)
X <- transposedXData[interX]
typeof(joinLip)
class(joinLip)
str(joinLip)
char <- as.character(joinLip)
as.character(char)
Y <- transposedYData[interY]
Y <- transposedYData[c("17351", "27432")]


cca.fit <- yacca::cca(X, Y)
cca.fit


helio.plot(cca.fit, x.name = "xLabel", y.name = "yLabel")




# PLS work. :)
library(pls)
yarn
oliveoil


model.matrix( ~ EIF2S3X + BAAT, data = transposedXData)
list(transposedXData, transposedYData)

colnames(transposedXData) <- paste("X", colnames(transposedXData), sep = ".")
colnames(transposedYData) <- paste("Y", colnames(transposedYData), sep = ".")

# WORKED! :D
library(pls)

Xmelt <- I(as.matrix(transposedXData))
class(Xmelt)
Ymelt <- I(as.matrix(transposedYData))

combined <- data.frame(X = I(Xmelt), Y = I(Ymelt))
combined <- data.frame(I(Xmelt),I(Ymelt))

PLSResultsFromMatrixInDF <- plsr(Y ~ X, data = combined )
summary(PLSResultsFromMatrixInDF)
plot(RMSEP(PLSResultsFromMatrixInDF), legendpos = "topright")
plot(PLSResultsFromMatrixInDF, ncomp = 10, asp = 1, line = TRUE)


combined <- data.frame(transposedXData, Y = I(transposedYData))



Xmelt <- I(as.matrix(transposedXData))
Ymelt <- I(as.matrix(transposedYData))
realData = list(Xmelt,Ymelt)
library(pls)
PLSResults <- plsr(Y.61205 + Y.61204 ~ X.EIF2S3X + X.ABCA1 + X.FABP6, data = combined )
summary(PLSResults)
plot(RMSEP(PLSResults), legendpos = "topright")

# Prediction plot.
plot(PLSResults, ncomp = 3, asp = 1, line = TRUE)
plot(PLSResults, plottype = "scores", comps = 1:3)
explvar(PLSResults)
validation="LOO"

plot(PLSResults, "loadings", comps = 1:2, legendpos = "topleft", labels = "numbers", xlab = "nm")
predict(PLSResults, ncomp = 2, newdata = transposedXData[c("X.TRA2B", "X.AATF")])

transposedXData[c("X.TRA2B", "X.AATF")]

pathToFileWithXData = "/home/koralgooll/doktorat/Rpackages/ONION/example/nm-transcriptomics.txt"
pathToFileWithYData = "/home/koralgooll/doktorat/Rpackages/ONION/example/nm-lipidomics-valid.txt"
XData <- readWithoutDuplicates(pathToFileWithXData)
YData <- readWithoutDuplicates(pathToFileWithYData)

transposedXData <- as.data.frame(t(XData))
transposedYData <- as.data.frame(t(YData))

nrow(transposedXData)


trainingRows <- ceiling(0.85 * nrow(transposedXData))
transposedXData[1:trainingRows,]
transposedXData[(trainingRows + 1):nrow(transposedXData),]


interX <- intersect(colnames(transposedXData), xNamesVector)
interY <- intersect(colnames(transposedYData), yNamesVector)
