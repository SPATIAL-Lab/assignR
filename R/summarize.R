jointP = function(pdR){
  
  if(!inherits(pdR, c("RasterStack", "RasterBrick", "SpatRaster"))){
    stop("input probability density map (pdR) should be a SpatRaster")
  }
  if(!inherits(pdR, "SpatRaster")){
    warning("raster objects are depreciated, transition to package terra")
    pdR = rast(pdR)
  }
  
  n = nlyr(pdR)
  result = pdR[[1]] * pdR[[2]]
  if(n > 2){
    for(i in seq_len(n)[-1:-2]){
      result = result * pdR[[i]]
    }
  }
  
  result = result / global(result, sum, na.rm = TRUE)[1, 1]
  names(result) = "Joint_Probability"
  p = options("scipen")
  on.exit(options(p))
  options(scipen = -2)
  plot(result)
  title("Joint Probability")
  return(result)
}

unionP = function(pdR){
  
  if(!inherits(pdR, c("RasterStack", "RasterBrick", "SpatRaster"))){
    stop("input probability density map (pdR) should be a SpatRaster")
  }
  if(!inherits(pdR, "SpatRaster")){
    warning("raster objects are depreciated, transition to package terra")
    pdR = rast(pdR)
  }
  
  result = (1 - pdR[[1]])
  n = nlyr(pdR)
  for(i in seq_len(n)[-1]){
    result = lapp(c(result, pdR[[i]]), fun = function(x, y){return(x*(1-y))})
  }
  
  plot(1-result)
  title("Union Probability")
  return(1-result)
}