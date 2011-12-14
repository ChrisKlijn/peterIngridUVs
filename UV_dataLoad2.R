# UV_dataLoad2.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Load IC50 data - second analysis pass
#              this analysis also contains experiments done in the
#              first analysis but 
# -------------------------------------------------------------------
# Run locally

# Libraries and functions
source('~/gitCodeChris/peterIngridUVs/UV_functions.R')

setwd("~/work/Matlab/Data/NKI Data/PeterIngridUV/rawData2/")
filesIC50 <- dir(pattern='.txt')
IC50List <- loadIC50files(filevect=filesIC50)

# Background subtraction and data shaping
IC50List.corr <- lapply(IC50List, correctIC50)

# Load sampleinfo

setwd("/home/klijn/data/smallproj/UVs")
sampleInfo <- read.delim('sampleinfo2.txt', stringsAsFactors=F, sep=',')
row.names(sampleInfo) <- sampleInfo$sampID

save(file='IC50DataRaw2.Rda', list=c('sampleInfo', 'IC50List.corr'))
