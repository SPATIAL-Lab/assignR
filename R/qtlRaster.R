qtlRaster = function(pdR, threshold, thresholdType = "area", 
                     genplot = TRUE, outDir = NULL){
  
  if(!inherits(pdR, c("RasterLayer", "RasterStack", "RasterBrick", "SpatRaster"))){
    stop("input probability density map (pdR) should be a SpatRaster")
  }
  if(!inherits(pdR, "SpatRaster")){
    warning("raster objects are depreciated, transition to package terra")
    pdR = rast(pdR)
  }
  if(!inherits(threshold, "numeric")){
    stop("threshold must be a number between 0 and 1 ")
  }
  if(length(threshold)[1] != 1){
    stop("threshold must be a number between 0 and 1 ")
  }
  if(threshold < 0 | threshold > 1){
    stop("threshold must be a number between 0 and 1")
  }
  if(thresholdType != "area" & thresholdType != "prob"){
    stop("thresholdType must be 'area' or 'prob'. See help page for 
         further information")
  }
  if(!inherits(genplot, "logical")) {
    message("genplot should be logical (T or F), using default = T")
    genplot = TRUE
  }
  if(!is.null(outDir)){
    if(!inherits(outDir, "character")){
      stop("outDir should be a character string")
    }
    if(!dir.exists(outDir)){
      message("outDir does not exist, creating")
      dir.create(outDir)
    }
  }
  
  result = pdR
  n = nlyr(result)  
  
  if(thresholdType == "prob"){
    if(threshold == 0){
      cut = rep(1, n)
    } else if(threshold == 1){
      cut = rep(0, n)
    } else{
      cut = double(n)
      pdR.values = values(pdR, na.rm = TRUE)
      for(i in seq_len(n)){
        pdR.value = sort(pdR.values[,i])
        k = length(pdR.value)
        left = 1
        right = k
        while((right-left) > 2){
          start = round(mean(c(left, right)))
          total = sum(pdR.value[start:k])
          if(total > threshold){
            left = start
          }
          if(total < threshold){
            right = start
          }
        }
        cut[i] = pdR.value[start]
      }
    }
    result = pdR > cut
    title1 = "probability"
  }
  
  if(thresholdType == "area"){
    if(threshold == 0){
      cut = rep(1, n)
    } else if(threshold == 1){
      cut = rep(0, n)
    } else{
      pdR.values = values(pdR, na.rm = TRUE)
      cut = apply(pdR.values, 2, quantile, probs = 1 - threshold)
    }
    result = pdR > cut
    title1 = "area"
  }
  
  names(result) = names(pdR)
  tls = character(n)
  if(n > 1){
    for(i in seq_len(n)){
      tls[i] = paste0("Top ", threshold*100, "% by ", title1, " for ", names(result)[i])
    }
  } else{
    tls = paste0("Top ", threshold*100, "% by ", title1)
  }
  
  if(genplot){
    for(i in seq_len(n)){
      plot(result[[i]], legend=FALSE)
      title(tls[i])
    }
  }
  if(!is.null(outDir)){
    pdf(paste0(outDir, "/qtlRaster_result.pdf"))
    for(i in seq_len(n)){
      plot(result[[i]], legend=FALSE)
      title(tls[i])
    }
    dev.off()
  }
  return(result)
}
