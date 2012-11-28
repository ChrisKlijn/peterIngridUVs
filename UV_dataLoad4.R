# UV_dataLoad2.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Load IC50 data - second analysis pass
#              this analysis also contains experiments done in the
#              first analysis but 
# -------------------------------------------------------------------
# Run on Macbook

# Libraries and functions
source('~/OldNKI/peterIngridUVs/UV_functions.R')

setwd('/Users/klijnc1/Projects/NKIProjects/UV/rawData')
filesIC50 <- dir(pattern='.txt')
IC50List <- loadIC50files(filevect=filesIC50)

# Background subtraction and data shaping
IC50List.corr <- lapply(IC50List, correctIC50)

# Load sampleinfo

setwd('/Users/klijnc1/Projects/NKIProjects/UV')
sampleInfo <- read.delim('sampleinfo4.txt', stringsAsFactors=F, sep='\t')
row.names(sampleInfo) <- sampleInfo$sampID

save(file='IC50DataRaw4.Rda', list=c('sampleInfo', 'IC50List.corr'))
