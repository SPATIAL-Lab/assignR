isoStack = function(..., clean = TRUE){

  r = list(...)
  
  if(class(r[[1]])[1] == "list"){
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
    if(!inherits(r[[i]], c("RasterBrick", "RasterStack"))){
      stop("each object in ... must be a RasterBrick or RasterStack")
    }
    if(nlayers(r[[i]]) != 2){
      stop("each isoscape must include two layers: mean and 1 sd")
    }
    if(is.na(proj4string(r[[i]]))) {
      stop("each isoscape must have valid coordinate reference system")
    }
  }
  
  if(compareRaster(r, rowcol = FALSE, res = TRUE, stopiffalse = FALSE)){
    #mask
    r = maskIso(r, n)
    
    #assign class
    class(r) = "isoStack"
    
    return(r)
  }
  
  if(clean == FALSE){
    stop("isoscape properties differ, clean set to FALSE")
  }
  
  #Get target proj
  c = crs(proj4string(r[[1]]))    

  if(!compareRaster(r, extent = FALSE, rowcol = FALSE, rotation = FALSE,
                    stopiffalse = FALSE)){
    for(i in 2:n){
      if(!compareCRS(r[[i]], c)){
        r[[i]] = projectRaster(r[[i]], crs = c)
      }
    }
  }
  
  #Get other target properties
  res.max = res(r[[1]])
  for(i in 2:n){
    res.max = pmin(res.max, res(r[[i]]))
  }
  
  e = extent(r[[1]])
  for(i in 2:n){
    e = intersect(e, extent(r[[i]]))
  }
  
  #Make raster target
  r.targ = raster(ext = e, crs = c, resolution = res.max)

  for(i in 1:n){
    if(!compareCRS(r[[i]], r.targ)){
      r[[i]] = projectRaster(r[[i]], r.targ)
    } else if(!compareRaster(r[[i]], r.targ, rowcol = FALSE, crs = FALSE,
                             res = TRUE, stopiffalse = FALSE)){
      r[[i]] = resample(r[[i]], r.targ)
    }
  }

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
    if(!inherits(i, c("RasterBrick", "RasterStack"))){
      stop("each object in r must be a RasterBrick or RasterStack")
    }
    if(nlayers(i) != 2){
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
