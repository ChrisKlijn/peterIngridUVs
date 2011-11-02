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

setwd("/home/klijn/data/smallproj/UVs/rawData")

load('IC50DataRaw.Rda')

fitList <- fitDualResponseCurve(targetFrame=IC50List.corr[[1]], controlFrame=IC50List.corr[[2]])

# png(file='testPlot.png', width=500, height=500)
plotDualResponseCurve(fitList)

# dev.off()

# Fit all curves and plot them

RMCE <- IC50List.corr[grep('RMCE', names(IC50List.corr))]
RMCE.fit <- lapply(RMCE, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
png(file="Figures/RMCEcurves.png", width=500, height=500)
  plot(RMCE.fit[[1]], type='n', main='RMCE dose response curves')
  lapply(RMCE.fit, plotSingleResponseCurve)
dev.off()

Brca1 <- IC50List.corr[grep('Brca1', names(IC50List.corr))]
Brca1.fit <- lapply(Brca1, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
png(file="Figures/Brca1curves.png", width=500, height=500)
plot(Brca1.fit[[1]], type='n', main='Brca1 dose response curves')
lapply(Brca1.fit, plotSingleResponseCurve)
dev.off()

UV <- IC50List.corr[grep('UV', names(IC50List.corr))]
UV.fit <- lapply(UV, function (x) { drm(response ~ conc, data=x, fct=LL.2())})
png(file="Figures/UVcurves.png", width=500, height=500)
plot(UV.fit[[1]], type='n', main='UV dose response curves')
lapply(UV.fit, plotSingleResponseCurve)
dev.off()

RMCE.IC50 <- lapply(RMCE.fit, function(x) {
  return(coef(x)['e:(Intercept)'])})
Brca1.IC50 <- lapply(Brca1.fit, function(x) {
  return(coef(x)['e:(Intercept)'])})
UV.IC50 <- lapply(UV.fit, function(x) {
  return(coef(x)['e:(Intercept)'])})
plot(RMCE.IC50, Brca1.IC50, pch=19, col=colors()[108], 
  main='IC50s of Brca1 and RMCE experiments', xlim=c(0,.3), ylim=c(0,.3))

allFit <- lapply(IC50List.corr, function (x) { drm(response ~ conc, 
  data=x, fct=LL.2())})

infoFrame <- getFitInfoFrame(allFit, sampleInfo)

qplot(data=subset(infoFrame, grepl('RMCE|Brca', sampID)), x=as.factor(series), y=IC50, color=type, geom='boxplot')

# normalize the IC50s to the RMCE

infoFrame.corr <- correctIC50(infoFrame)

# Some plots
qplot(data=infoFrame.corr, x=as.factor(series), y=IC50corr, color=type, geom='boxplot')
qplot(data=infoFrame.corr, x=IC50corr, color=type, geom='boxplot')
qplot(data=infoFrame.corr, x=type, y=IC50corr, color=type, size=RSE, geom='jitter')

# Check the poor fits
poorFits <- subset(infoFrame.corr, RSE > .10)
pdf('Figures/poorFits.pdf', width=5, height=5)
for (i in 1:nrow(poorFits)) {
  plotResponseCurve(allFit[poorFits$sampID][[i]], title=poorFits$sampID[i])
}
dev.off()

# Pathogenicity score?

path <- c(mu=mean(subset(infoFrame.corr, pathogenic==1)$IC50corr),
  sigma=sd(subset(infoFrame.corr, pathogenic==1)$IC50corr))
nonpath <- c(mu=mean(subset(infoFrame.corr, pathogenic==0)$IC50corr),
  sigma=sd(subset(infoFrame.corr, pathogenic==0)$IC50corr))
infoFrame.corr$zPath <- (infoFrame.corr$IC50corr - path['mu'])/
  path['sigma']
infoFrame.corr$zNonPath <- (infoFrame.corr$IC50corr - nonpath['mu'])/
  nonpath['sigma']



infoFrame.corr$pPath <- pnorm(q=infoFrame.corr$distPath + path['mu'], 
  mean=path['mu'], sd=path['sigma'], lower.tail=F)
infoFrame.corr$pNonPath <- pnorm(q=infoFrame.corr$distNonPath + 
  nonpath['mu'], mean=nonpath['mu'], sd=nonpath['sigma'], lower.tail=F)


##


a <- rnorm(n=10, .5)
b <- rnorm(n=10, 2.2)

infoFrame$type <- ordered(infoFrame$type, c('RMCE', 'hBrca1', 'UV'))
normList <- split(infoFrame, infoFrame$series)
infoFrame.corrLin <- 
  Reduce(f=rbind, lapply(normList, normalizeIC50linear))

tempCorrIC50 <- function(Data) {

  subData <- subset(Data, grepl('RMCE|Brca', sampID))
  tempLM <- lm(data=subData, IC50 ~ type)
  subData$IC50 <- (subData$IC50 - coef(tempLM)[1]) / coef(tempLM)[2]
  
  return(subData)

}

testCorr <- Reduce(f=rbind, lapply(testList, tempCorrIC50))

qplot(data=testCorr, x=as.factor(series), y=IC50, color=type, geom='boxplot')
plot(0,0, type='n', xlim=c(-1,1), ylim=c(-1,1))
