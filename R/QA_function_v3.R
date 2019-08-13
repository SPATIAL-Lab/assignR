QA <- function(isoscape, known, valiStation, valiTime, setSeed = T){

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
    stop("known must have same coordinate reference system as isoscape")
  } else if(ncol(known@data) != 1){
    stop("known must include a 1-column data frame containing only the isotope values")
  }

  if(!valiStation<nrow(known)){
    stop("valiStation must be smaller than the number of known-origin stations in known")
  }
  if(setSeed == T){
    set.seed(100)
  }
  
  rowLength <- nrow(known)
  val_stations <- sort(sample(1:rowLength,valiStation,replace = F))
  for (i in 1:(valiTime-1)){
    val_stations <- rbind(val_stations,sort(sample(1:rowLength,valiStation,replace = F)))
  }

  stationNum4model <- rowLength - valiStation
  prption_byProb <- matrix(0, valiTime, 99) # accuracy by checking top percentage by cumulative prob.
  prption_byArea <- matrix(0, valiTime, 99) # accuracy by checking top percentage by area
  pd_bird_val <- matrix(0, valiTime, valiStation) # pd value for each validation location
  precision <- list() # precision
  
  # create progress bar
  pb <- txtProgressBar(min = 0, max = valiTime, style = 3)
  
  for (i in 1:valiTime){
    bird_val <- known[val_stations[i,],]
    bird_model <- known[-val_stations[i,],]
    rescale <- isOrigin::calRaster(bird_model, isoscape, sdMethod = 1, genplot = F, savePDF = F, verboseLM = F)
    pd <- isOrigin::pdRaster(rescale, unknown = data.frame(row.names(bird_val@data), bird_val@data[,1]), genplot = F, saveFile = F)

    # pd value for each validation location
    for(m in 1:nlayers(pd)){
      pd_bird_val[i, m] <- raster::extract(pd[[m]], bird_val[m,], method = "bilinear")
    }

    xx <- seq(1, 99, 1) ## 1 to 99

    # total area
    Tarea <- length(na.omit(pd[[1]][]))

    # accuracy and precision by checking top percentage by cumulative prob.
    precision[[i]] <- matrix(0, 99, valiStation) # precision
    for(j in xx){
      qtl <- isOrigin::qtlRaster(pd, threshold = j/100, pdf = F, thresholdType = 1,genplot = F)
      prption_byProb[i, j] <- 0
      for(k in 1:nlayers(qtl)){
        prption_byProb[i, j] <- prption_byProb[i, j] +
          raster::extract(qtl[[k]], bird_val[k,], method = "bilinear")
        precision[[i]][j, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
      }
    }

    for(k in 1:nlayers(qtl)){
      precision[[i]][j, k] <- sum(na.omit(qtl[[k]][]))/Tarea # precision
    }

    # accuracy by checking top percentage by cumulative area
    for(n in xx){
      qtl <- isOrigin::qtlRaster(pd, threshold = n/100, pdf = F, thresholdType = 2,genplot = F)
      prption_byArea[i, n] <- 0
      for(k in 1:nlayers(qtl)){
        prption_byArea[i, n] <- prption_byArea[i, n] +
          raster::extract(qtl[[k]], bird_val[k,], method = "bilinear")
      }
    }
   
    #update progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }

  random_prob_density=1/length(na.omit(getValues(isoscape[[1]])))

  result <- list(val_stations, pd_bird_val, prption_byArea, prption_byProb, precision, random_prob_density)
  names(result) <- c("val_stations", "pd_bird_val", "prption_byArea", "prption_byProb", "precision", "random_prob_density")
  class(result) <- "QA"
  return(result)
}
