# UV_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Sandbox code
# -------------------------------------------------------------------

# Medoid

# Libraries and functions
library(reshape)
library(drc)
source('~/codeChris/smallProjects/peterIngridUV/UV_functions.R')

# Load raw files

setwd("/home/klijn/data/smallproj/UVs/rawData")
filesIC50 <- dir(pattern='.txt')
IC50List <- loadIC50files(filevect=filesIC50)

# Background subtraction and data shaping
IC50List.corr <- lapply(IC50List, correctIC50)

# Load sampleinfo

setwd("/home/klijn/data/smallproj/UVs")
sampleInfo <- read.delim('sampleinfo.csv', stringsAsFactors=F, sep=',')


fitList <- fitDualResponseCurve(targetFrame=IC50List.corr[[1]], controlFrame=IC50List.corr[[2]])

# png(file='testPlot.png', width=500, height=500)
plotDualResponseCurve(fitList)

# dev.off()

# test control hybs

RMCE <- IC50List.corr[grep('RMCE', names(IC50List.corr))]
RMCE.fit <- lapply(RMCE, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
x11()
plot(RMCE.fit[[1]], type='n', main='RMCE dose response curves')
lapply(RMCE.fit, plotSingleResponseCurve)

# test control hybs

Brca1 <- IC50List.corr[grep('Brca1', names(IC50List.corr))]
Brca1.fit <- lapply(Brca1, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
x11()
plot(Brca1.fit[[1]], type='n', main='Brca1 dose response curves')
lapply(Brca1.fit, plotSingleResponseCurve)

RMCE.IC50 <- lapply(RMCE.fit, function(x) {
  return(coef(x)['e:(Intercept)'])})
Brca1.IC50 <- lapply(Brca1.fit, function(x) {
  return(coef(x)['e:(Intercept)'])})
plot(RMCE.IC50, Brca1.IC50, pch=19, col=colors()[108], 
  main='IC50s of Brca1 and RMCE experiments', xlim=c(0,.3), ylim=c(0,.3))

  