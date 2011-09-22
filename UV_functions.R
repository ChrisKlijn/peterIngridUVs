# UV_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: BRCA1 Unclassified Variants
# Description: Functions to support BRCA1 UV analyses
# -------------------------------------------------------------------

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

fitResponseCurve <- function(targetFrame, controlFrame, fitFunc=LL.2()) {
  
  targetFit <- drm(response ~ conc, data=targetFrame, fct=fitFunc)
  controlFit <- drm(response ~ conc, data=controlFrame, fct=fitFunc)

  return(list(target=targetFit, 
    control=controlFit, 
    targetIC50=coef(targetFit)['e:(Intercept)'],
    controlIC50=coef(controlFit)['e:(Intercept)']))

}

plotResponseCurve <- function (fitList) {
  
  plot(fitList$target, col='red', pch=19, type='average', cex=.8)
  plot(fitList$target, col='red', pch=19, type='bars', add=T, lty='blank')
  abline(v=fitList$targetIC50, col='red')
  plot(fitList$control, col='blue', pch=19, type='average', cex=.8, add=T)
  plot(fitList$control, col='blue', pch=19, type='bars', add=T, lty='blank')
  abline(v=fitList$controlIC50, col='blue')

}
