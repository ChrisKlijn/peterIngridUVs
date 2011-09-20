# UV_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Sandbox code
# -------------------------------------------------------------------

setwd('~/Desktop/IngridTest')

loadIC50files <- function (filevect) {
  
  fileList <- vector(mode='list', length=length(filevect))
  names(fileList) <- filevect

  for (i in filevect) {
    
    fileList[[i]] <- read.delim(file=i, stringsAsFactors=F)
          
  }

  names(fileList) <- gsub('.txt', '', names(fileList))
  return(fileList)
}

correctIC50 <- function(IC50Rawframe) {
  
  require(reshape)

  IC50Rawframe.corr <- IC50Rawframe[,2:ncol(IC50Rawframe)] - 
    mean(IC50Rawframe[,'no.cells'])
  IC50Rawframe.corr <- IC50Rawframe.corr / mean(IC50Rawframe.corr[,1])
  colnames(IC50Rawframe.corr) <- gsub('X', '', colnames(IC50Rawframe.corr))

  tempFrame <- melt(IC50Rawframe.corr)
  colnames(tempFrame) <- c('conc', 'response')
  tempFrame$conc <- as.numeric(levels(tempFrame$conc))[tempFrame$conc]

  return(tempFrame)  
    
}

filesIC50 <- dir(pattern='.txt')

IC50List <- loadIC50files(filevect=filesIC50)
IC50List.corr <- lapply(IC50List, correctIC50)

a <- drm(response ~ conc, data=IC50List.corr[[1]], fct=LL.2())
b <- drm(response ~ conc, data=IC50List.corr[[2]], fct=LL.2())

png(file='testPlot.png', width=500, height=500)
plot(a, col='blue')
plot(b, col='red', add=T)
abline(v=coef(a)[2], col='blue')
abline(v=coef(b)[2], col='red')
legend('bottomleft', legend=names(IC50List.corr), fill=c('red','blue'))
dev.off()

