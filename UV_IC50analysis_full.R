# UV_IC50analysis_full.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Analysis of IC50 data using all measurements
# -------------------------------------------------------------------

# Run locally

# Libraries and functions
library(reshape)
library(drc)
library(ggplot2)
source('~/gitCodeChris/peterIngridUV/UV_functions.R')

setwd("~/work/Matlab/Data/NKI Data/PeterIngridUV/")

load('IC50DataRaw.Rda')

# Fit the data

allFit <- lapply(IC50List.corr, function (x) { drm(response ~ conc, 
  data=x, fct=LL.2())})

infoFrame <- getFitInfoFrame(allFit, sampleInfo)

# normalize the IC50s to the RMCE

# Simple RMCE scaling
infoFrame.corr <- normalizeIC50(infoFrame)

# Linear model for RMCE to Brca1
infoFrame$type <- ordered(infoFrame$type, c('RMCE', 'hBrca1', 'UV'))
normList <- split(infoFrame, infoFrame$series)
infoFrame.corrLin <- 
  Reduce(f=rbind, lapply(normList, normalizeIC50linear))

# Fit the data, excluding the first measurement

IC50List.corr.ex1 <- lapply(IC50List.corr, function (x) {
  return(x[x$conc != 0,])})
allFit.ex1 <- lapply(IC50List.corr.ex1, function (x) { 
  drm(response ~ conc, data=x, fct=LL.2())})
infoFrame.ex1 <- getFitInfoFrame(allFit.ex1, sampleInfo)
# normalize the IC50s to the RMCE
infoFrame.ex1.corr <- normalizeIC50(infoFrame.ex1)

# Pvalues 

infoFrame.corrLin <- calcBayesProb(infoFrame.corrLin)

# Export to text

write.table(x=infoFrame.corrLin, file="IC50analysis.txt", col.names=TRUE,  quote=FALSE, row.names=FALSE, sep='\t')

# Plots

png(file="Figures/full_RMCE_Brca_boxplot.png", width=800, height=500)
qplot(data=subset(infoFrame, grepl('RMCE|Brca', sampID)), x=as.factor(series), y=IC50, color=type, geom='boxplot')
dev.off()

# Some more plots
png(file="Figures/full_RMCE_Brca_UV_norm_boxplot.png", 
  width=800, height=500)
qplot(data=infoFrame.corr, x=as.factor(series), y=IC50corr, color=type, geom='boxplot')
dev.off()

png(file="Figures/full_density_corrIC50.png", width=500, height=500)
qplot(data=infoFrame.corr, x=IC50corr, fill=type, alpha=I(.5), 
  geom='density')
dev.off()

png(file="Figures/full_jitter_corrIC50.png", width=500, height=500)
qplot(data=infoFrame.corr, x=type, y=IC50corr, color=type, size=RSE, geom='jitter')
dev.off()

## Plots for the linear normalization

png(file="Figures/full_RMCE_Brca_UV_normLin_boxplot.png", 
  width=800, height=500)
qplot(data=infoFrame.corrLin, x=as.factor(series), y=IC50, color=type, geom='boxplot')
dev.off()

png(file="Figures/full_density_corrLinIC50.png", width=500, height=500)
qplot(data=infoFrame.corrLin, x=IC50, fill=type, alpha=I(.5), 
  geom='density')
dev.off()

png(file="Figures/full_jitter_corrLinIC50.png", width=500, height=500)
qplot(data=infoFrame.corrLin, x=type, y=IC50, color=type, size=RSE, geom='jitter')
dev.off()

png(file="Figures/full_pvaljit_corrLinIC50.png", width=500, height=500)
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


# Density of the RSE

