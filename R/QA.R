QA = function(known, isoscape, bySite = TRUE, 
              valiStation = 1, valiTime = 50, by = 2, 
              mask = NULL, setSeed = TRUE, name = NULL){

  #check bySite
  if(!is.logical(bySite)){
    stop("bySite must be logical")
  }

  #check that known is valid and has defined, correct CRS
  col_site = NULL
  if(!(class(known)[1] %in% c("subOrigData", "SpatialPointsDataFrame"))) {
    stop("known must be a subOrigData or SpatialPointsDataFrame object")
  }
  if(class(known)[1] == "subOrigData"){
    if(is.null(known$data) || is.null(known$marker)){
      stop("missing information in subOrigData known")
    }
    col_m = match(known$marker, names(known$data))
    col_sd = match(paste0(known$marker, ".sd"), names(known$data))
    if(is.na(col_m) | is.na(col_sd)){
      stop("cannot match marker to data table in subOrigData known")
    }
    if(bySite){
      col_site = match("Site_ID", names(known$data))
      if(is.na(col_site)){
        stop("no Site_ID field in known; provide Site_IDs or use bySite = FALSE")
      }
    }
    known = known$data
    known@data = known@data[, c(col_m, col_sd, col_site)]
  }else{
    if(ncol(known@data) < 2){
      if(ncol(known@data) == 1){
        warning("use of known with 1 data column is depreciated; known 
              should include a data frame containing the measured 
              isotope values (col 1) and 1 sd uncertainty (col 2); 
              assuming equal uncertainty for all samples")
        known@data[,2] = rep(0.0001, nrow(known@data))
      } else{
        stop("known must include a data frame containing the measured 
              isotope values (col 1) and 1 sd uncertainty (col 2)")
      }
    }
    if(any(!is.numeric(known@data[,1]), !is.numeric(known@data[,2]))){
      stop("known must include a data frame containing the measured 
           isotope values (col 1) and 1 sd uncertainty (col 2)")
    }
    warning("user-provided known; assuming measured isotope value and 1 sd
            uncertainty are contained in columns 1 and 2, respectively")
    if(bySite){
      col_site = match("Site_ID", names(known@data))
      if(is.na(col_site)){
        stop("no Site_ID field in known; provide Site_IDs or use bySite = FALSE")
      }
      known@data = known@data[, c(1, 2, col_site)]
    }
  }
  if(any(is.na(known@data[, 1])) || 
     any(is.nan(known@data[, 1])) || 
     any(is.null(known@data[, 1]))){
    stop("Missing values detected in known values")
  }
  if(any(is.na(known@data[, 2])) || 
     any(is.nan(known@data[, 2])) || 
     any(is.null(known@data[, 2]))){
    stop("Missing values detected in known uncertainties")
  }
  if(any(known@data[, 2] == 0)){
    stop("zero values found in known uncertainties")
  }
  if(!is.null(col_site)){
    if(any(is.na(known@data[, col_site])) || 
       any(is.nan(known@data[, col_site])) || 
       any(is.null(known@data[, col_site]))){
      stop("Missing values detected in sites field")
    }
  }
  if(is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } 
  if(proj4string(known) != proj4string(isoscape)){
    known = spTransform(known, crs(isoscape))
    warning("known was reprojected")
  } 
  if(nrow(known) < 10){
    warning("there are fewer than 10 known samples")
  }
  if(nrow(known) < 3){
    stop("QA requires at least 3 known samples")
  }
  
  #check isoscape
  if (class(isoscape)[1] == "RasterStack" | 
      class(isoscape)[1] == "RasterBrick") {
    if (is.na(sp::proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
  } else {
    stop("isoscape should be a RasterStack or RasterBrick")
  }
  
  #check valiStation
  if(bySite){
    if(valiStation > (length(unique(known$Site_ID)) - 3)){
      stop("for bySite = TRUE valiStation must be 3 or more smaller 
        than the number of unique sites in known")
    }
  } else{
    if(valiStation > (nrow(known) - 3)){
      stop("for bySite = FALSE valiStation must be 3 or more smaller 
        than the number of samples in known")
    }
  }
  
  #check valiTime
  if(valiTime < 2){
    stop("valiTime must be an integer greater than 1")
  }
  
  #check by
  if(!(as.integer(by) == by) || by < 1 || by > 25){
    stop("by must be an integer between 1 and 25")
  }
  
  #check mask
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
  
  #check name
  if(!is.null(name)){
    if(class(name)[1] != "character"){
      stop("name must be a character string")
    }
  }
  
  #check setSeed
  if(!is.logical(setSeed)){
    stop("setSeed must be logical")
  }
  
  if(setSeed){
    set.seed(100)
  }
  
  if(bySite){
    rowLength = length(unique(known$Site_ID))
    ids = known$Site_ID
  } else{
    rowLength = nrow(known)
    ids = seq_len(rowLength)
  }
  val_stations = sort(sample(ids, valiStation, replace = FALSE))
  for (i in seq_len(valiTime)[-1]){
    val_stations = rbind(val_stations, 
                          sort(sample(ids, valiStation, 
                                      replace = FALSE)))
  }

  xx = seq(1, 101, by)
  if(tail(xx, 1) != 101){
    xx = c(xx, 101)
  }
  prption_byProb = matrix(0, valiTime, length(xx))   
  prption_byArea = matrix(0, valiTime, length(xx)) 
  pd_v = matrix(0, valiTime, valiStation) 
  precision = list() 
  
  # create progress bar
  pb = txtProgressBar(min = 0, max = valiTime, style = 3)
  
  for (i in seq_len(valiTime)){
    if(bySite){
      v = known[known$Site_ID %in% val_stations[i,],]
      m = known[-(known$Site_ID %in% val_stations[i,]),]
    } else{
      v = known[val_stations[i,],]
      m = known[-val_stations[i,],]
    }
    
    class(m) = "QAData"
    rescale = calRaster(m, isoscape, mask, genplot = FALSE, 
                                  verboseLM = FALSE)
    
    pd = pdRaster(rescale, 
                  unknown = data.frame(row.names(v@data), v@data[,1]), 
                  genplot = FALSE)

    # pd value for each validation sample or site
    pd_temp = double(nlayers(pd))
    for(j in seq_len(nlayers(pd))){
      pd_temp[j] = extract(pd[[j]], v[j,])
    }
    if(bySite){
      for(j in seq(valiStation)){
        pd_v[i, j] = mean(pd_temp[v$Site_ID == val_stations[i, j]])
      }
    } else{
      pd_v[i,] = pd_temp
    }

    # total area
    Tarea = length(na.omit(pd[[1]][]))

    # spatial precision and accuracy by checking top percentage by cumulative prob.
    precision[[i]] = matrix(0, length(xx), valiStation)
    for(j in seq_along(xx)){
      qtl = qtlRaster(pd, threshold = (xx[j]-1)/100, 
                                thresholdType = "prob", 
                               genplot = FALSE)
      rv_temp = double(nlayers(qtl))
      pre_temp = double(nlayers(qtl))
      for(k in seq_len(nlayers(qtl))){
        rv_temp[k] = extract(qtl[[k]], v[k,])
        pre_temp[k] = sum(na.omit(qtl[[k]][]))/Tarea
      }
      if(bySite){
        rv_sm = double(valiStation)
        for(k in seq(valiStation)){
          rv_sm[k] = mean(rv_temp[v$Site_ID == val_stations[i, k]], 
                          na.rm = TRUE)
          precision[[i]][j, k] = mean(pre_temp[v$Site_ID == 
                                                 val_stations[i, k]])
        }
        if(any(!is.na(rv_sm))){
          prption_byProb[i, j] = mean(rv_sm, na.rm = TRUE)
        } else{
          prption_byProb[i, j] = NA
        }
      } else{
        if(any(!is.na(rv_sm))){
          prption_byProb[i, j] = sum(rv_temp, na.rm = TRUE) /
            sum(!is.na(rv_temp))
        } else{
          prption_byProb[i, j] = NA
        }
        precision[[i]][j,] = pre_temp
      }
    }

    # sensitivity by checking top percentage by cumulative area
    for(j in seq_along(xx)){
      qtl = qtlRaster(pd, threshold = (xx[j]-1)/100, 
                                thresholdType = "area", 
                               genplot = FALSE)
      rv_temp = double(nlayers(qtl))
      for(k in seq_len(nlayers(qtl))){
        rv_temp[k] = extract(qtl[[k]], v[k,])
      }
      if(bySite){
        rv_sm = double(valiStation)
        for(k in seq(valiStation)){
          rv_sm[k] = mean(rv_temp[v$Site_ID == val_stations[i, k]], 
                          na.rm = TRUE)
        }
        if(any(!is.na(rv_sm))){
          prption_byArea[i, j] = mean(rv_sm, na.rm = TRUE)
        } else{
          prption_byArea[i, j] = NA
        }
      } else{
        if(any(!is.na(rv_sm))){
          prption_byArea[i, j] = sum(rv_temp, na.rm = TRUE) /
            sum(!is.na(rv_temp))
        } else{
          prption_byArea[i, j] = NA
        }
      }
    }
   
    #update progress bar
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }

  random_prob_density=1/length(na.omit(getValues(isoscape[[1]])))

  result = list(name, val_stations, pd_v, prption_byArea, 
                prption_byProb, precision, random_prob_density, by)
  names(result) = c("name", "val_stations", "pd_val", 
                    "prption_byArea", "prption_byProb", "precision",
                    "random_prob_density", "by")
  class(result) = "QA"
  return(result)
}
