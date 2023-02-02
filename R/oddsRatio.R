oddsRatio = function(pdR, inputP){
  
  if(!inherits(pdR, c("RasterLayer", "RasterStack", "RasterBrick", "SpatRaster"))){
    stop("input probability density map (pdR) should be a SpatRaster")
  }
  if(!inherits(pdR, "SpatRaster")){
    warning("raster objects are depreciated, transition to package terra")
    pdR = rast(pdR)
  }

  if(inherits(inputP, "SpatialPoints")){
    if(is.na(proj4string(inputP))){
      stop("inputP must have coord. ref.")
    }
    if(proj4string(inputP) != crs(pdR, proj = TRUE)){
      inputP = spTransform(inputP, crs(pdR))
      message("inputP was reprojected")
    }
    
    n = length(inputP)
    extrVals = extract(pdR, vect(inputP))[,-1]
    if(any(is.na(extrVals))){
      stop("one or more points have probability NA")
    }
    gmax = matrix(rep(global(pdR, "max", na.rm = TRUE)[,1], n), 
                  nrow = n, byrow = TRUE)
    gmin = matrix(rep(global(pdR, "min", na.rm = TRUE)[,1], n), 
                  nrow = n, byrow = TRUE)
    result2 = data.frame(ratioToMax = extrVals/gmax, 
                         ratioToMin = extrVals/gmin)
    if(n == 1){
      result = result2
    }
    else if(n == 2){
      if(nlyr(pdR) > 1){
        result1 = (extrVals[1,]/(1-extrVals[1,])) / 
          (extrVals[2,]/(1-extrVals[2,]))
      } else {
        result1 = (extrVals[1]/(1-extrVals[1])) / (extrVals[2]/(1-extrVals[2]))
      }
      result = list(oddsRatio = result1, ratioToMaxMin = result2)
      names(result) = c("P1/P2 odds ratio", "Odds relative to the max/min pixel")
      row.names(result[[2]]) = c("P1", "P2")
    }
    else{
      stop("input points (inputP) should be one or two points")
    }
  }
  
  if(inherits(inputP, "SpatialPolygons")){
    if(length(inputP) != 2){
      stop("input polygons (inputP) should be two polygons")
    }
    if(is.na(proj4string(inputP))){
      stop("inputP must have coord. ref.")
    }
    if(proj4string(inputP) != crs(pdR, proj = TRUE)){
      inputP = spTransform(inputP, crs(pdR))
      message("inputP was reprojected")
    }
    
    extrVals = extract(pdR, vect(inputP), "sum", na.rm = TRUE)
    if(any(extrVals[, -1] == 0)){
      stop("No values in P1 and/or P2")
    }
    
    result1 = (extrVals[1, -1]/(1-extrVals[1, -1])) /
      (extrVals[2, -1]/(1-extrVals[2, -1]))
    result2 = length(cells(crop(pdR, vect(inputP[1,])))) / 
      length(cells(crop(pdR, vect(inputP[2,]))))
    result = list(oddsRatio = result1, polygonCellRatio = result2)
    names(result) = c("P1/P2 odds ratio", "Ratio of numbers of cells in two polygons")
  }

  return(result)
}
