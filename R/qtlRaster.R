qtlRaster <- function(pdR, threshold, thresholdType = "area", genplot = TRUE, savePDF = FALSE){
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
  if(thresholdType != "area" & thresholdType != "prob"){
    stop("thresholdType must be 'area' or 'prob'. See help page for further information")
  }
  if(class(genplot) != "logical"){
    stop("genplot must be logical (T/F)")
  }
  if(class(savePDF) != "logical"){
    stop("pdf must be logical (T/F)")
  }
  
  result <- pdR
  n = raster::nlayers(result)  
  if(thresholdType == "prob"){
    for(i in 1:n){
      if(threshold == 0){
        cut = 1
      } else if(threshold == 1){
        cut = 0
      } else{
        pdR.values <- stats::na.omit(raster::getValues(pdR[[i]]))
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
        cut = pdR.values[start]        
      }
      if(n == 1){
        result <- pdR[[i]] > cut
      }else{
        result[[i]] <- pdR[[i]] > cut
      }
    }
    title1 <- "probability"
  }
  
  if(thresholdType == "area"){
    for(i in 1:n){
      if(threshold == 0){
        cut = 1
      } else if(threshold == 1){
        cut = 0
      } else{
        pdR.values <- stats::na.omit(raster::getValues(pdR[[i]]))
        k <- length(pdR.values)
        cut <- sort(pdR.values)[round((1-threshold)*k)]
      }
      if(n == 1){
        result <- pdR[[i]] > cut
      }else{
        result[[i]] <- pdR[[i]]>cut
      }
    }
    title1 <- "area"
  }
  
  names(result) <- names(pdR)
  tls = character(n)
  if(n > 1){
    for(i in 1:n){
      tls[i] = paste0("Top ", threshold*100, "% quantile by ", title1, " for ", names(result)[i])
    }
  } else{
    tls = paste0("Top ", threshold*100, "% quantile by ", title1)
  }
  
  if(genplot){
    for(i in 1:n){
      graphics::plot(result[[i]], legend=FALSE)
      graphics::title(tls[i])
    }
  }
  if(savePDF){
    grDevices::pdf("qtlRaster_result.pdf")
    for(i in 1:n){
      graphics::plot(result[[i]], legend=FALSE)
      graphics::title(tls[i])
    }
    grDevices::dev.off()
  }
  return(result)
}