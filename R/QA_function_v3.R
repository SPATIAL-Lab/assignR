QA <- function(isoscape, known, valiStation, valiTime, mask = NULL, setSeed = TRUE){

  #check that isoscape is valid and has defined CRS
  if (class(isoscape) == "RasterStack" | class(isoscape) == "RasterBrick") {
    if (is.na(proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
  } else {
    stop("isoscape should be a RasterStack or RastrBrick")
  }

  #check that known is valid and has defined, correct CRS
  if (class(known) != "SpatialPointsDataFrame") {
    stop("known should be a SpatialPointsDataFrame, see help page of calRaster function")
  } else if (is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } else if(proj4string(known) != proj4string(isoscape)){
    known = spTransform(known, crs(isoscape))
    warning("known was reprojected")
  } else if(ncol(known@data) != 1){
    stop("known must include a 1-column data frame containing only the isotope values")
  }

  if(!valiStation<nrow(known)){
    stop("valiStation must be smaller than the number of known-origin stations in known")
  }
  
  if(!is.null(mask)) {
    if(class(mask) == "SpatialPolygonsDataFrame" || class(mask) == "SpatialPolygons"){
      if(is.na(proj4string(mask))){
        stop("mask must have coordinate reference system")
      }
      if(proj4string(mask) != proj4string(isoscape)){
        mask <- spTransform(mask, crs(isoscape))
        warning("mask was reprojected")
      }
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }
  
  if(setSeed == TRUE){
    set.seed(100, sample.kind = "Rejection")
  }
  
  rowLength <- nrow(known)
  val_stations <- sort(sample(1:rowLength,valiStation,replace = FALSE))
  for (i in 1:(valiTime-1)){
    val_stations <- rbind(val_stations,sort(sample(1:rowLength,valiStation,replace = FALSE)))
  }

  stationNum4model <- rowLength - valiStation
  prption_byProb <- matrix(0, valiTime, 99) # accuracy by checking top percentage by cumulative prob.
  prption_byArea <- matrix(0, valiTime, 99) # accuracy by checking top percentage by area
  pd_v <- matrix(0, valiTime, valiStation) # pd value for each validation location
  precision <- list() # precision
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = valiTime, style = 3)
  
  for (i in 1:valiTime){
    v <- known[val_stations[i,],]
    m <- known[-val_stations[i,],]
    rescale <- assignR::calRaster(m, isoscape, mask, genplot = FALSE, verboseLM = FALSE)
    pd <- assignR::pdRaster(rescale, unknown = data.frame(row.names(v@data), v@data[,1]), genplot = FALSE)

    # pd value for each validation location
    for(j in 1:nlayers(pd)){
      pd_v[i, j] <- raster::extract(pd[[j]], v[j,], method = "bilinear")
    }

    xx <- seq(1, 99, 1) ## 1 to 99

    # total area
    Tarea <- length(na.omit(pd[[1]][]))

    # spatial precision and accuracy by checking top percentage by cumulative prob.
    precision[[i]] <- matrix(0, 99, valiStation) # precision
    for(j in xx){
      qtl <- assignR::qtlRaster(pd, threshold = j/100, savePDF = FALSE, thresholdType = 1, genplot = FALSE)
      prption_byProb[i, j] <- 0
      for(k in 1:nlayers(qtl)){
        rv = raster::extract(qtl[[k]], v[k,], method = "bilinear")
        if(!is.na(rv)){
          prption_byProb[i, j] <- prption_byProb[i, j] + rv
        }
        precision[[i]][j, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
      }
    }

    for(k in 1:nlayers(qtl)){
      precision[[i]][j, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
    }

    # sensitivity by checking top percentage by cumulative area
    for(n in xx){
      qtl <- assignR::qtlRaster(pd, threshold = n/100, savePDF = FALSE, thresholdType = 2, genplot = FALSE)
      prption_byArea[i, n] <- 0
      for(k in 1:nlayers(qtl)){
        rv = raster::extract(qtl[[k]], v[k,], method = "bilinear")
        if(!is.na(rv)){
          prption_byArea[i, n] <- prption_byArea[i, n] + rv
        }
      }
    }
   
    #update progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }

  random_prob_density=1/length(na.omit(getValues(isoscape[[1]])))

  result <- list(val_stations, pd_v, prption_byArea, prption_byProb, precision, random_prob_density)
  names(result) <- c("val_stations", "pd_val", "prption_byArea", "prption_byProb", "precision", "random_prob_density")
  class(result) <- "QA"
  return(result)
}
