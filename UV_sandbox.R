# UV_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Sandbox code
# -------------------------------------------------------------------

# WD - local for now

setwd('~/Desktop/IngridTest')

# Libraries and functions
library(reshape)
library(drc)
source('~/codeChris/smallProjects/peterIngridUV/UV_functions.R')

filesIC50 <- dir(pattern='.txt')

IC50List <- loadIC50files(filevect=filesIC50)
IC50List.corr <- lapply(IC50List, correctIC50)

fitList <- fitResponseCurve(targetFrame=IC50List.corr[[1]], controlFrame=IC50List.corr[[2]])

# png(file='testPlot.png', width=500, height=500)
plotResponseCurve(fitList)

# dev.off()
