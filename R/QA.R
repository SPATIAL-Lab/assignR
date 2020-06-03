QA = function(isoscape, known, 
              valiStation = ceiling(length(known)*0.1), 
               valiTime = 50, by = 2, mask = NULL, 
              setSeed = TRUE, name = NULL){

  #check that isoscape is valid and has defined CRS
  if (class(isoscape)[1] == "RasterStack" | 
      class(isoscape)[1] == "RasterBrick") {
    if (is.na(sp::proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
  } else {
    stop("isoscape should be a RasterStack or RasterBrick")
  }

  #check that known is valid and has defined, correct CRS
  if (class(known)[1] != "SpatialPointsDataFrame") {
    stop("known should be a SpatialPointsDataFrame, 
         see help page of calRaster function")
  }
  if(any(is.na(known@data[,1])) || any(is.nan(known@data[,1])) || 
     any(is.null(known@data[,1]))){
    stop("Missing values detected in known")
  }
  if (is.na(sp::proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } 
  if(sp::proj4string(known) != sp::proj4string(isoscape)){
    known = sp::spTransform(known, raster::crs(isoscape))
    warning("known was reprojected")
  } else if(ncol(known@data) != 1){
    stop("known must include a 1-column data frame containing only the 
         isotope values")
  }
  if(nrow(known) < 10){
    warning("there are fewer than 10 known samples")
  }
  if(nrow(known) < 3){
    stop("QA requires at least 3 known samples")
  }
  if(!valiStation<nrow(known)){
    stop("valiStation must be smaller than the number of known-origin 
         stations in known")
  }
  if(valiTime < 2){
    stop("valiTime must be an integer greater than 1")
  }
  if(!(as.integer(by) == by) || by < 1 || by > 25){
    stop("by must be an integer between 1 and 25")
  }
  if(!is.null(mask)) {
    if(class(mask)[1] == "SpatialPolygonsDataFrame" || 
       class(mask)[1] == "SpatialPolygons"){
      if(is.na(sp::proj4string(mask))){
        stop("mask must have coordinate reference system")
      }
      if(sp::proj4string(mask) != sp::proj4string(isoscape)){
        mask = sp::spTransform(mask, raster::crs(isoscape))
        warning("mask was reprojected")
      }
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }
  if(!is.null(name)){
    if(class(name)[1] != "character"){
      stop("name must be a character string")
    }
  }
  
  if(setSeed == TRUE){
    set.seed(100)
  }
  
  rowLength = nrow(known)
  val_stations = sort(sample(seq_len(rowLength), valiStation, replace = FALSE))
  for (i in seq_len(valiTime)[-1]){
    val_stations = rbind(val_stations, 
                          sort(sample(seq_len(rowLength), valiStation, 
                                      replace = FALSE)))
  }

  stationNum4model = rowLength - valiStation
  xx = seq(1, 101, by)
  if(tail(xx, 1) != 101){
    xx = c(xx, 101)
  }
  prption_byProb = matrix(0, valiTime, length(xx))   
  prption_byArea = matrix(0, valiTime, length(xx)) 
  pd_v = matrix(0, valiTime, valiStation) 
  precision = list() 
  
  # create progress bar
  pb = utils::txtProgressBar(min = 0, max = valiTime, style = 3)
  
  for (i in seq_len(valiTime)){
    v = known[val_stations[i,],]
    m = known[-val_stations[i,],]
    rescale = assignR::calRaster(m, isoscape, mask, genplot = FALSE, 
                                  verboseLM = FALSE)
    pd = assignR::pdRaster(rescale, 
                            unknown = data.frame(row.names(v@data), 
                                                 v@data[,1]), genplot = FALSE)

    # pd value for each validation location
    for(j in seq_len(raster::nlayers(pd))){
      pd_v[i, j] = raster::extract(pd[[j]], v[j,])
    }

    # total area
    Tarea = length(stats::na.omit(pd[[1]][]))

    # spatial precision and accuracy by checking top percentage by cumulative prob.
    precision[[i]] = matrix(0, length(xx), valiStation) # precision
    for(j in seq_along(xx)){
      qtl = assignR::qtlRaster(pd, threshold = (xx[j]-1)/100, 
                                thresholdType = "prob", 
                               genplot = FALSE)
      prption_byProb[i, j] = 0
      for(k in seq_len(raster::nlayers(qtl))){
        rv = raster::extract(qtl[[k]], v[k,])
        if(!is.na(rv)){
          prption_byProb[i, j] = prption_byProb[i, j] + rv
        }
        precision[[i]][j, k] = sum(stats::na.omit(qtl[[k]][]))/Tarea 
      }
    }

    # sensitivity by checking top percentage by cumulative area
    for(n in seq_along(xx)){
      qtl = assignR::qtlRaster(pd, threshold = (xx[n]-1)/100, 
                                thresholdType = "area", 
                               genplot = FALSE)
      prption_byArea[i, n] = 0
      for(k in seq_len(raster::nlayers(qtl))){
        rv = raster::extract(qtl[[k]], v[k,])
        if(!is.na(rv)){
          prption_byArea[i, n] = prption_byArea[i, n] + rv
        }
      }
    }
   
    #update progress bar
    Sys.sleep(0.1)
    utils::setTxtProgressBar(pb, i)
  }

  random_prob_density=1/length(stats::na.omit(raster::getValues(isoscape[[1]])))

  result = list(name, val_stations, pd_v, prption_byArea, 
                prption_byProb, precision, random_prob_density, by)
  names(result) = c("name", "val_stations", "pd_val", 
                    "prption_byArea", "prption_byProb", "precision",
                    "random_prob_density", "by")
  class(result) = "QA"
  return(result)
}
