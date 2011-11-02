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

fitDualResponseCurve <- function(targetFrame, controlFrame, 
  fitFunc=LL.2()) {
  
  targetFit <- drm(response ~ conc, data=targetFrame, fct=fitFunc)
  controlFit <- drm(response ~ conc, data=controlFrame, fct=fitFunc)

  return(list(target=targetFit, 
    control=controlFit, 
    targetIC50=coef(targetFit)['e:(Intercept)'],
    controlIC50=coef(controlFit)['e:(Intercept)']))

}

plotResponseCurve <- function(fitObj, title='test') {
   
   plot(fitObj, col='red', pch=19, type='obs', cex=.8, 
    main=title)
   plot(fitObj, col='red', pch=19, type='bars', add=T)

}

plotDualResponseCurve <- function (fitList) {
  
  plot(fitList$target, col='red', pch=19, type='average', cex=.8)
  plot(fitList$target, col='red', pch=19, type='bars', add=T, lty='blank')
  abline(v=fitList$targetIC50, col='red')
  plot(fitList$control, col='blue', pch=19, type='average', cex=.8, add=T)
  plot(fitList$control, col='blue', pch=19, type='bars', add=T, lty='blank')
  abline(v=fitList$controlIC50, col='blue')

}

plotSingleResponseCurve <- function(fitObj) {
  
  plot(fitObj, col='red', pch=19, type='average', cex=.8, add=T)
  plot(fitObj, col='red', pch=19, type='bars', add=T, lty='blank')
  abline(v=fitObj$targetIC50['e:(Intercept)'], col='red')
}

getFitInfoFrame <- function (fitList, sampleInfo) {
  
  all.IC50 <- unlist(lapply(fitList, function(x) {
    return(coef(x)['e:(Intercept)'])}))
  names(all.IC50) <- gsub('[.]e[:][(]Intercept[)]', '', names(all.IC50))
  all.RSE <- unlist(lapply(fitList, function(x) {
    return(summary(x)$rseMat[1,1])}))
  
  all.IC50.frame <- data.frame(IC50=all.IC50, RSE=all.RSE)
  rownames(all.IC50.frame) <- names(all.IC50)

  returnFrame <- merge(all.IC50.frame, sampleInfo, by='row.names')
  returnFrame$type <- gsub('[_].*$|[0-9]{1,2}[.].*$', '',
    returnFrame$sampID)

  return(returnFrame)
}

normalizeIC50 <- function(IC50Frame) {
  
  .corrSeriesFrame <- function(seriesFrame) {
    
    RMCEval <- mean(seriesFrame$IC50[
      grepl('RMCE', seriesFrame$sampID)])
    seriesFrame$IC50corr <- seriesFrame$IC50 / RMCEval

    return(seriesFrame)
  }

  # Split the data according to series

  seriesList <- split(IC50Frame, f=IC50Frame$series)
  seriesList.corr <- lapply(seriesList, .corrSeriesFrame)

  return(Reduce(f=rbind, x=seriesList.corr))

}

normalizeIC50linear <- function(IC50Frame) {

  subData <- subset(IC50Frame, grepl('RMCE|Brca', sampID))
  tempLM <- lm(data=subData, IC50 ~ type)
  IC50Frame$IC50 <- (IC50Frame$IC50 - coef(tempLM)[1]) / 
    coef(tempLM)[2]
  
  return(IC50Frame)
  
}