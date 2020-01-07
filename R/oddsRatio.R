oddsRatio <- function(pdR, inputP){
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }

  if(class(inputP) == "SpatialPoints" || class(inputP) == "SpatialPointsDataFrame"){
    if(is.na(sp::proj4string(inputP))){
      stop("inputP must have coord. ref.")
    }
    if(sp::proj4string(inputP) != sp::proj4string(pdR)){
      inputP = sp::spTransform(inputP, raster::crs(pdR))
      warning("inputP was reprojected")
    }
    
    n <- length(inputP)
    extrVals <- raster::extract(pdR, inputP)
    result2 <- data.frame(ratioToMax = extrVals/raster::maxValue(pdR), ratioToMin = extrVals/raster::minValue(pdR))
    if(n == 1){
      result = result2
    }
    else if(n == 2){
      if(class(pdR) == "RasterStack" | class(pdR) == "RasterBrick"){
        result1 <- (extrVals[1,]/(1-extrVals[1,])) / (extrVals[2,]/(1-extrVals[2,]))
      } else {
        result1 <- (extrVals[1]/(1-extrVals[1])) / (extrVals[2]/(1-extrVals[2]))
      }
      result <- list(oddsRatio = result1, ratioToMaxMin = result2)
      names(result) <- c("P1/P2 odds ratio", "Odds relative to the max/min pixel")
      row.names(result[[2]]) = c("P1", "P2")
    }
    else{
      stop("input points (inputP) should be one or two points")
    }
  }
  
  if(class(inputP) == "SpatialPolygons" || class(inputP) == "SpatialPolygonsDataFrame"){
    if(length(inputP) != 2){
      stop("input polygons (inputP) should be two polygons")
    }
    if(is.na(sp::proj4string(inputP))){
      stop("inputP must have coord. ref.")
    }
    if(sp::proj4string(inputP) != sp::proj4string(pdR)){
      inputP = sp::spTransform(inputP, raster::crs(pdR))
      warning("inputP was reprojected")
    }
    
    extrVals <- raster::extract(pdR, inputP)
    if(nlayers(pdR) > 1){
      extrVals.p1 <- colSums(extrVals[[1]])
      extrVals.p2 <- colSums(extrVals[[2]])
    } else {
      extrVals.p1 = sum(extrVals[[1]])
      extrVals.p2 = sum(extrVals[[2]])
    }
    result1 <- (extrVals.p1/(1-extrVals.p1))/(extrVals.p2/(1-extrVals.p2))
    result2 <- raster::ncell(raster::crop(pdR, inputP[1,]))/raster::ncell(raster::crop(pdR, inputP[2,]))
    result <- list(oddsRatio = result1, polygonCellRatio = result2)
    names(result) <- c("P1/P2 odds ratio", "Ratio of numbers of cells in two polygons")
  }

  return(result)
}
