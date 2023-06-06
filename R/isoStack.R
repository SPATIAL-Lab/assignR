isoStack = function(..., clean = TRUE){

  r = list(...)
  
  if(inherits(r[[1]], "list")){
    r = unlist(r, recursive = FALSE)
  }
  
  if((!inherits(r, "list")) | length(r) < 2){
    stop("... must be a list containing multiple isoscapes")
  }
  n = length(r)
  
  for(i in 1:n){
    if(inherits(r[[i]], "rescale")){
      r[[i]] = r[[i]]$isoscape.rescale
    }
    if(!inherits(r[[i]], c("RasterBrick", "RasterStack", "SpatRaster"))){
      stop("each object in ... must be a SpatRaster")
    }
    if(inherits(r[[i]], c("RasterBrick", "RasterStack"))){
      r[[i]] = rast(r[[i]])
      warning("raster objects are depreciated, transition to package terra")
    }
    if(nlyr(r[[i]]) != 2){
      stop("each isoscape must include two layers: mean and 1 sd")
    }
    if(crs(r[[i]]) == "") {
      stop("each isoscape must have valid coordinate reference system")
    }
  }

  #projections
  for(i in 2:n){
    if(!same.crs(r[[i]], r[[1]])){
      if(clean){
        r[[i]] = project(r[[i]], crs(r[[1]]))
      } else{
        stop("isoscape projections differ, clean set to FALSE")
      }
    }
  }  
  
  #check other properties
  res.flag = FALSE
  ext.flag = FALSE
  res.max = r[[1]]@ptr$res
  ext.min = ext(r[[1]])
  for(i in 2:n){
    if(!identical(r[[i]]@ptr$res, r[[1]]@ptr$res)) res.flag = TRUE
    if(!identical(ext(r[[i]]), ext(r[[1]]))) ext.flag = TRUE
    
    res.max = pmin(res.max, r[[i]]@ptr$res)
    ext.min = intersect(ext.min, ext(r[[i]]))
  }
  
  #fix other properties
  if(res.flag | ext.flag){
    if(clean){
      #Make raster target
      r.targ = rast(ext = ext.min, resolution = res.max,
                    crs = crs(r[[1]]))
      
      for(i in 1:n){
        if(!compareGeom(r[[i]], r.targ, rowcol = FALSE, crs = FALSE,
                                 res = TRUE, stopOnError = FALSE)){
          r[[i]] = resample(r[[i]], r.targ)
        }
      }
    } else{
      stop("isoscape properties differ, clean set to FALSE")
    }
  }

  #common mask
  r = maskIso(r, n)
  
  #assign class
  class(r) = "isoStack"
  
  return(r)
}

plot.isoStack = function(x, ...){
  
  if(!inherits(x, "isoStack")){
    stop("plot.isoStack needs isoStack object")
  }
  
  if(length(x) < 2){
    stop("isoStack must include at least 2 isoscapes")
  }
  
  for(i in x){
    if(!inherits(i, c("RasterBrick", "RasterStack", "SpatRaster"))){
      stop("each object in r must be a SpatRaster")
    }
    if(nlyr(i) != 2){
      stop("each isoscape must include two layers: mean and 1 sd")
    }
  }
  
  for(i in x){
    plot(i)
  }
}

maskIso = function(r, n){
  #Create mask
  m = r[[1]]
  for(i in 2:n){
    m = m * r[[i]]
  }
  m = m[[1]] * m[[2]]
  
  #Apply mask
  for(i in 1:n){
    r[[i]] = mask(r[[i]], m)
  }
  
  return(r)
}
