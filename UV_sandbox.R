# UV_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Sandbox code
# -------------------------------------------------------------------

# WD - medoid

setwd("/home/klijn/data/smallproj/UVs")

# Libraries and functions
library(reshape)
library(drc)
source('~/codeChris/smallProjects/peterIngridUV/UV_functions.R')

# Load files

setwd("/home/klijn/data/smallproj/UVs/rawData")
filesIC50 <- dir(pattern='.txt')

IC50List <- loadIC50files(filevect=filesIC50)
IC50List.corr <- lapply(IC50List, correctIC50)

fitList <- fitDualResponseCurve(targetFrame=IC50List.corr[[1]], controlFrame=IC50List.corr[[2]])

# png(file='testPlot.png', width=500, height=500)
plotDualResponseCurve(fitList)

# dev.off()

# test control hybs

a <- IC50List.corr[grep('RMCE', names(IC50List.corr))]
b <- lapply(a, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
x11()
plot(b[[1]], type='n')
lapply(b, plotSingleResponseCurve)
