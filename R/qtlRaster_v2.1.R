qtlRaster <- function(pdR, threshold, thresholdType, genplot = T, pdf = F){
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }
  if(class(threshold) != "numeric"){
    stop("threshold must be one number between 0 and 1 ")
  }
  if(length(threshold) != 1){
    stop("threshold must be one number between 0 and 1 ")
  }
  if(threshold < 0 | threshold > 1){
    stop("threshold must be one number between 0 and 1")
  }
  if(length(thresholdType) != 1){
    stop("thresholdType must be 1 or 2. See help page for further information")
  }
  if(thresholdType != 1 & thresholdType != 2){
    stop("thresholdType must be 1 or 2. See help page for further information")
  }
  if(class(genplot) != "logical"){
    stop("genplot must be logical (T/F)")
  }
  if(class(pdf) != "logical"){
    stop("pdf must be logical (T/F)")
  }
  result <- pdR
  if(thresholdType == 1){
    for(i in 1:nlayers(pdR)){
      pdR.values <- na.omit(getValues(pdR[[i]]))
      pdR.values <- sort(pdR.values)
      k <- length(pdR.values)
      left <- 1
      right <-  k
      while((right-left) > 2){
        start <- round((left+right)/2)
        total <- sum(pdR.values[start:k])
        if(total > threshold){
          left <- start
        }
        if(total < threshold){
          right <- start
        }
      }
      if(nlayers(pdR) == 1){
        result <- pdR[[i]] > pdR.values[start]
      }else{
        result[[i]] <- pdR[[i]] > pdR.values[start]
      }
    }
    title1 <- "by Cumulative Probability"
  }
  if(thresholdType == 2){
    for(i in 1:nlayers(pdR)){
      pdR.values <- na.omit(getValues(pdR[[i]]))
      k <- length(pdR.values)
      cut <- sort(pdR.values)[round((1-threshold)*k)]
      if(nlayers(pdR) == 1){
        result <- pdR[[i]] > cut
      }else{
        result[[i]] <- pdR[[i]]>cut
      }
    }
    title1 <- "by Area"
  }
  names(result) <- names(pdR)
  if(genplot){
    for(i in 1:nlayers(result)){
      plot(result[[i]])
      title(paste0("Top ", threshold*100, "% ", title1, " for ", names(result)[i]))
    }
  }
  if(pdf){
    pdf("qtlRaster_result.pdf")
    for(i in 1:nlayers(result)){
      plot(result[[i]])
      title(paste0("Top ", threshold*100, "% ", title1, " for ", names(result)[i]))
    }
    dev.off()
  }
  return(result)
}
