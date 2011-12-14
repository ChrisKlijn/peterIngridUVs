# UV_IC50analysis_2.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Analysis of IC50 data using all measurements.
#              This is the second analysis done on 8 series of data
# -------------------------------------------------------------------

# Run locally

# Libraries and functions
library(reshape)
library(drc)
library(ggplot2)
source('~/gitCodeChris/peterIngridUVs/UV_functions.R')

setwd("~/work/Matlab/Data/NKI Data/PeterIngridUV/")

load('IC50DataRaw2.Rda')

# Fit the data

allFit <- lapply(IC50List.corr, function (x) { drm(response ~ conc, 
  data=x, fct=LL.2())})
infoFrame <- getFitInfoFrame(allFit, sampleInfo)

# normalize the IC50s to the RMCE

# Linear model for RMCE to Brca1
infoFrame$type <- ordered(infoFrame$type, c('RMCE', 'hBrca1', 'UV'))
normList <- split(infoFrame, infoFrame$series)
infoFrame.corrLin <- 
  Reduce(f=rbind, lapply(normList, normalizeIC50linear))

# Pvalues 
infoFrame.corrLin <- calcBayesProb(infoFrame.corrLin)

# Export to text
write.table(x=infoFrame.corrLin, file="IC50analysis2.txt", col.names=TRUE,  quote=FALSE, row.names=FALSE, sep='\t')

# --------------- Plotting -------------------

## Plots to check the normalization

png(file="Figures/2_RMCE_Brca_UV_normLin_boxplot.png", 
  width=800, height=500)
qplot(data=infoFrame.corrLin, x=as.factor(series), y=IC50, color=type, geom='boxplot')
dev.off()

png(file="Figures/2_density_corrLinIC50.png", width=500, height=500)
qplot(data=infoFrame.corrLin, x=IC50, fill=type, alpha=I(.5), 
  geom='density')
dev.off()

png(file="Figures/2_jitter_corrLinIC50.png", width=500, height=500)
qplot(data=infoFrame.corrLin, x=type, y=IC50, color=type, size=RSE, geom='jitter')
dev.off()

png(file="Figures/2_pvaljit_corrLinIC50.png", width=500, height=500)
qplot(data=infoFrame.corrLin, x=type, y=IC50corr, color=type, shape=Pval > .05, size=RSE, geom='jitter')
dev.off()

# Comparison lineair to non-normalized and simple normalized

png(file="Figures/norm_comparison.png", width=800, height=500)
par(mfrow=c(1,3))
boxplot(infoFrame.corrLin$IC50 ~ infoFrame.corrLin$type, 
  main='normal IC50')
boxplot(infoFrame.corr$IC50corr ~ infoFrame.corrLin$type, 
  main='Simple norm')
boxplot(infoFrame.corrLin$IC50corr ~ infoFrame.corrLin$type, 
  main='Lineair norm')
dev.off()

# Check the poor fits
poorFits <- subset(infoFrame.corr, RSE > .10)
pdf('Figures/poorFits.pdf', width=5, height=5)
for (i in 1:nrow(poorFits)) {
  plotResponseCurve(allFit[poorFits$sampID][[i]], title=poorFits$sampID[i])
}
dev.off()

# Compare RSE for inlcusion or exlcusion of the first measurement

pdf('Figures/comparison_firstMeasurement.pdf', width=5, height=5)
qplot(data=infoFrame.corr, x=RSE, fill=type, alpha=I(.5), geom='histogram')
qplot(data=infoFrame.ex1.corr, x=RSE, fill=type, alpha=I(.5), geom='histogram')
dev.off()

# Check the poor fits
poorFits.ex1 <- subset(infoFrame.ex1.corr, RSE > .10)
pdf('Figures/poorFits_ex1.pdf', width=5, height=5)
for (i in 1:nrow(poorFits.ex1)) {
  plotResponseCurve(allFit.ex1[poorFits.ex1$sampID][[i]], title=poorFits.ex1$sampID[i])
}
dev.off()
