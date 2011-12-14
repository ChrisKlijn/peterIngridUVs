# UV_dataLoad.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Load IC50 data
# -------------------------------------------------------------------

# Run locally

# Libraries and functions
library(reshape)
source('~/gitCodeChris/peterIngridUV/UV_functions.R')

# Load raw files

setwd("~/work/Matlab/Data/NKI Data/PeterIngridUV/rawData/")
filesIC50 <- dir(pattern='.txt')
IC50List <- loadIC50files(filevect=filesIC50)

# Background subtraction and data shaping
IC50List.corr <- lapply(IC50List, correctIC50)

# Load sampleinfo

setwd("/home/klijn/data/smallproj/UVs")
sampleInfo <- read.delim('sampleinfo.csv', stringsAsFactors=F, sep=',')
row.names(sampleInfo) <- sampleInfo$sampID

save(file='IC50DataRaw.Rda', list=c('sampleInfo', 'IC50List.corr'))
