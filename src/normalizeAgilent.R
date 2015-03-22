# ===============================================================================================
# ========================= NORMALIZE DATA; AGILENT ========================================
# ===============================================================================================
options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)
#require{agilp}
# R 3.0.2

# package for normalisation agilent data:
#source("http://bioconductor.org/biocLite.R")
#biocLite("agilp") 
library(agilp)

# ------------- function: AAProcess ----------------------------
# Extract and save the median expression data from Agilent scanner arrays,
# averaging replicate probes. The raw (unprocessed) expression data is frst extracted 
# from the scanner fles which are in tab delimited .txt format.

inputdir<-args[1] 
#outputdir<-("G:/PROJEKTY2012/CYFRONET/Ekspertyza/R_obliczenia/Normalizacja_agilent/agilent_raw/out")
outputdir<-args[2]
raw<-file.path(outputdir,"raw","", fsep = .Platform$file.sep)
AAProcess(input = inputdir, output = raw, s = 9)

# ------------- function: filenamex ----------------------------
# Copy extracted raw file names to template. This example makes a list of files in the folder 
# agilp/extdata/raw and saves it in a file called names.txt (tab delimited) in the folder agilp/output.
#filenamex(input=raw,output=outputdir)

# ------------- function: Baseline ----------------------------
# The AALoess normalisation function (see below) requires a `baseline' # file which contains the mean value for each probe from a number of dierent 
# arrays against which to normalise each new array. Baseline generates such 
# a baseline from a set of arrays.

#baselines
template<-file.path(outputdir,"names.txt", fsep = .Platform$file.sep)
template
inputdir
# Baseline for K(-) 15d
outputbase<-file.path(outputdir,"base", "baseline_K(-) 15d.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=6,B=9,input=raw, baseout=outputbase, t = template)

# Baseline for K(+) 15d
outputbase<-file.path(outputdir,"base", "baseline_K(+) 15d.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=2,B=5,input=raw, baseout=outputbase, t = template)

# Baseline for K(-) 30d
outputbase<-file.path(outputdir,"base", "baseline_K(-) 30d.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=14,B=17,input=raw, baseout=outputbase, t = template)

# Baseline for K(+) 30d
outputbase<-file.path(outputdir,"base", "baseline_K(+) 30d.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=10,B=13,input=raw, baseout=outputbase, t = template)

# Baseline for K(-) 48h
outputbase<-file.path(outputdir,"base", "baseline_K(-) 48h.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=22,B=25,input=raw, baseout=outputbase, t = template)

# Baseline for K(+) 48h
outputbase<-file.path(outputdir,"base", "baseline_K(+) 48h.txt", fsep = .Platform$file.sep)
Baseline(NORM="LOG",allfiles="TRUE",r=2,A=18,B=21,input=raw, baseout=outputbase, t = template)

# ------------- function: AALoess ----------------------------
# Normalises a set of gene expression data fles using LOESS

output<-file.path(outputdir,"output", "", fsep = .Platform$file.sep)
output
# AALoess for K(+) 48h
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(+) 48h.txt", fsep = .Platform$file.sep)
baselinedir
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE") 

# AALoess for K(-) 48h
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(-) 48h.txt", fsep = .Platform$file.sep)
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE") 

# AALoess for K(+) 15d
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(+) 15d.txt", fsep = .Platform$file.sep)
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE") 

# AALoess for K(-) 15d
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(-) 15d.txt", fsep = .Platform$file.sep)
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE") 

# AALoess for K(+) 30d
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(+) 30d.txt", fsep = .Platform$file.sep)
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE") 

# AALoess for K(-) 30d
#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
baselinedir<-file.path(outputdir,"base","baseline_K(-) 30d.txt", fsep = .Platform$file.sep)
AALoess(input=raw, output=output, baseline = baselinedir, LOG="TRUE")
#unlink(paste(le.path(system.le(package=\agilp"),\output",\"),\*.*",sep=\"), recursive=FALSE) #To remove these fles again and empty the output directory 


# ------------- function: IDswop ----------------------------
# Maps expression data across different bioinformatic identifiers.
# It is used files produced by biomaRt,
mapped<-file.path(outputdir,"map", "", fsep = .Platform$file.sep)
annotation<-args[3]
IDswop(input=output,output=mapped,annotation=annotation,source_ID="ProbeID",target_ID="ensembl_peptide_id", ERR=1)


# ------------- function: Equaliser  ----------------------------
# Takes files of raw data selects all entries common to all 
# input files can be output of the IDswop

#inputdir<-file.path(system.file(package="agilp"),"extdata","raw","", fsep = .Platform$file.sep)
Equaliser(input = mapped, output = output)


