wDist = function(pdR, sites){
  
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }
  if(is.na(proj4string(pdR))){
    stop("pdR must have coord. ref.")
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
  
  #make space
  ns = length(sites)
  db = data.frame(character(ns))
  for(i in 1:16){
    db = cbind(db, numeric(ns))
  }
  d.dens = list()
  b.dens = list()
  n = names(pdR)
  
  p = function(y, w){
    Position(function(x) x >= y, w)
  }
  
  for(i in 1:ns){
    pdSP = rasterToPoints(pdR[[i]], spatial = TRUE)
    
    d = distGeo(pdSP, sites[i,])
    b = bearing(sites[i,], pdSP)
    w = pdSP@data[,1]
    d.dens[[i]] = density(d, weights = w)
    b.dens[[i]] = density(b, weights = w)    
    
    #record weighted mean of distance distribution
    db[i, 1] = n[i]
    db[i, 2] = weighted.mean(d, w)
    
    #find and record quantiles within weighted distance distribution
    dw = cbind(d, w)
    dw = dw[order(d),]
    for(j in 2:nrow(dw)){
      dw[j, 2] = dw[j, 2] + dw[j-1, 2]
    }

    qts = sapply(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), p, w = dw[, 2])
    db[i, 3:9] = dw[qts, 1]
    
    bw = cbind(b, w)

    #find minimum weight value in bearing data to establish 'break'
    bbins = seq(-180, 170, by = 10)
    mp = 1
    for(j in bbins){
      pbin = sum(bw[bw[,1] >= j & bw[,1] < j + 10, 2])
      if(pbin < mp){
        mp = pbin
        mpb = j
      }
    }
    
    #re-reference bearing data to break
    bw[, 1] = bw[, 1] - mpb
    for(j in seq_along(bw[,1])){
      if(bw[j, 1] < 0){
        bw[j, 1] = bw[j, 1] + 360
      }
    }
    
    #weighted mean, re-referenced
    db[i, 10] = weighted.mean(bw[, 1], bw[, 2])
    
    #find and record quantiles within weighted bearing distribution
    bw = bw[order(bw[, 1]),]
    for(j in 2:nrow(bw)){
      bw[j, 2] = bw[j, 2] + bw[j-1, 2]
    }
    qts = sapply(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), p, w = bw[, 2])
    db[i, 11:17] = bw[qts, 1]
    
    #rectify bearings
    db[i, 10:17] = db[i, 10:17] + mpb
    for(j in 10:17){
      if(db[i, j] >= 180){
        db[i, j] = db[i, j] - 360
      }
    }
  }
  
  names(db) = c("sampleID", "wMeanDist", "w05Dist", "w10Dist", "w25Dist",
                "w50Dist", "w75Dist", "w90Dist", "w95Dist", "wMeanBear",
                "w05Bear", "w10Bear", "w25Bear", "w50Bear", "w75Bear",
                "w90Bear", "w95Bear")
  names(d.dens) = names(b.dens) = n
  
  wd = list(stats = db, d.dens = d.dens, b.dens = b.dens)
  class(wd) = "wDist"
  return(wd)
}
