wDist = function(pdR, sites){
  
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }
  if(proj4string(pdR) != "+proj=longlat +ellps=WGS84 +no_defs"){
    pdR = projectRaster(pdR, crs = "+proj=longlat +ellps=WGS84 +no_defs")
  }
  
  if(class(sites)[1] == "SpatialPoints" || 
     class(sites)[1] == "SpatialPointsDataFrame"){
    if(length(sites) != nlayers(pdR)){
      stop("sites and pdR have different lenghts; wDist requires one site per pdR layer")
    }
    if(is.na(proj4string(sites))){
      stop("sites must have coord. ref.")
    }
    if(proj4string(sites) != proj4string(pdR)){
      sites = spTransform(sites, crs(pdR))
    }
  } else{
    stop("sites should be a SpatialPoints object")
  }
  
  
  for(i in seq_along(sites)){
    pdSP = rasterToPoints(pdR[[i]])
    xy = coordi
  }
}