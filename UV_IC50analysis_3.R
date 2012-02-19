# UV_IC50analysis_2.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Analysis of IC50 data using all measurements.
#              This is the third analysis done, resulting in the
#              analysis that will end up in the paper.
#              Locations have been changed to perform the analyis on
#              my new computer.
# -------------------------------------------------------------------

# Run locally

# Libraries and functions
library(reshape)
library(drc)
library(ggplot2)
source('~/OldNKI/peterIngridUVs/UV_functions.R')

setwd('/Users/klijnc1/Projects/NKIProjects/UV')

load('IC50DataRaw3.Rda')

# Remove excluded experiments from sampleInfo and data

IC50List.corr[sampleInfo$sampID[sampleInfo$exclude == 1]] <- NULL
sampleInfo <- sampleInfo[!(sampleInfo$exclude == 1),]

# Check if the info and data for all UVs is present

sum(sampleInfo$sampID %in% names(IC50List.corr)) == nrow(sampleInfo)
sum(names(IC50List.corr) %in% sampleInfo$sampID) == 
  length(IC50List.corr)

# Fit the data

allFit <- lapply(IC50List.corr, function (x) { drm(response ~ conc, 
  data=x, fct=LL.2())})
infoFrame <- getFitInfoFrame(allFit, sampleInfo)

# Plot fits of all data

pdf(file='Figures/allFits_run3.pdf', width=7, height=7)
for (i in 1:length(allFit)) {
  plotResponseCurve(allFit[[i]], names(allFit)[i])
  tempInfo <- infoFrame[infoFrame$sampID == names(allFit)[i],]
  text(x=.8, y=1, labels=
    paste('RSE =', signif(tempInfo$RSE, digits=4)))
}
dev.off()

# normalize the IC50s to the RMCE

# Linear model for RMCE to Brca1
infoFrame$type <- ordered(infoFrame$type, c('RMCE', 'hBrca1', 'UV'))
normList <- split(infoFrame, infoFrame$series)
infoFrame.corrLin <- 
  Reduce(f=rbind, lapply(normList, normalizeIC50linear))

# Pvalues 
infoFrame.corrLin <- calcBayesProb(infoFrame.corrLin)
infoFrame.corrLin$PvalCorr <- p.adjust(infoFrame.corrLin$Pval, 
  method='BH')
infoFrame.corrLin$PvalPath <- infoFrame.corrLin$PvalCorr > 0.05

# Adjust the pathogenic score for plotting
infoFrame.corrLin$pathogenic[is.na(infoFrame.corrLin$pathogenic)] <- 2
infoFrame.corrLin$pathogenic <- 
  relevel(x=as.factor(infoFrame.corrLin$pathogenic), ref='1')

# Some checks

with(infoFrame.corrLin, table(status, series))
with(subset(infoFrame.corrLin, status == 'bad_duplo'), 
  table(UVID, PvalPath))
with(subset(infoFrame.corrLin, status == 'bad_fit'), 
  table(UVID, PvalPath))

# Export to text
write.table(x=infoFrame.corrLin, file="IC50analysis3.txt", col.names=TRUE,  quote=FALSE, row.names=FALSE, sep='\t')

# --------------- Plotting -------------------

## Plots to check the normalization

# Check the RSE per experiment
qplot(data=infoFrame.corrLin, x=as.factor(series), y=RSE, color=type, geom='boxplot')
ggsave("Figures/3_RSEcheck.png")
ggsave("Figures/3_RSEcheck.eps")

#-------------------------
# Paper plots all fits

qplot(data=infoFrame.corrLin, x=as.factor(series), 
  y=IC50corr, color=type, geom='boxplot') + 
  theme_bw() +
  opts(title='All 212 Experiments')
ggsave("Figures/3_IC50corr_check_all.png")
ggsave("Figures/3_IC50corr_check_all.pdf")

qplot(data=infoFrame.corrLin, x=type, y=IC50corr, 
  color=as.factor(pathogenic), shape=Pval > .05, size=I(3), 
  position=position_jitter(w=.3, h=0)) +
  theme_bw() +
  opts(title='All 212 Experiments')
ggsave("Figures/3_jitterplotIC50_all.png")
ggsave("Figures/3_jitterplotIC50_all.pdf")

qplot(data=infoFrame.corrLin, x=IC50corr, fill=type, alpha=I(.5), 
  geom='density') + theme_bw() +
  opts(title='All 212 Experiments')
ggsave("Figures/3_densityplotIC50_all.png")
ggsave("Figures/3_densityplotIC50_all.pdf")

#-------------------------
# Paper plots exclude = 0 fits

infoFrame.corrLinTopFit <- 
  infoFrame.corrLin[infoFrame.corrLin$exclude == 0,]

qplot(data=infoFrame.corrLinTopFit, x=as.factor(series), 
  y=IC50corr, color=type, geom='boxplot') + 
  theme_bw() +
  opts(title='Only include = 0 Fits')
ggsave("Figures/3_IC50corr_check_topfit.png")
ggsave("Figures/3_IC50corr_check_topfit.pdf")

qplot(data=infoFrame.corrLinTopFit, x=type, y=IC50corr, 
  color=as.factor(pathogenic), shape=Pval > .05, size=I(3), 
  position=position_jitter(w=.3, h=0)) +
  theme_bw() +
  opts(title='Only include = 0 Fits')
ggsave("Figures/3_jitterplotIC50_all_topfit.png")
ggsave("Figures/3_jitterplotIC50_all_topfit.pdf")

qplot(data=infoFrame.corrLinTopFit, x=IC50corr, fill=type, alpha=I(.5), 
  geom='density') + theme_bw() +
  opts(title='Only include = 0 Fits')
ggsave("Figures/3_densityplotIC50_topfit.png")
ggsave("Figures/3_densityplotIC50_topfit.pdf")
