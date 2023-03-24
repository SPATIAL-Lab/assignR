QA = function(known, isoscape, bySite = TRUE, valiStation = 1, 
              valiTime = 50, recal = TRUE, by = 2, prior = NULL, 
              mask = NULL, setSeed = TRUE, name = NULL){

  #space to handle messages and warnings
  mstack = wstack = character(0)
  addm = function(cnd){
    mstack <<- append(mstack, cnd$message)
    cnd_muffle(cnd)
  }
  addw = function(cnd){
    wstack <<- append(wstack, cnd$message)
    cnd_muffle(cnd)
  }
  
  #check bySite
  if(!is.logical(bySite)){
    stop("bySite must be logical")
  }

  #check isoscape and set ni number of isotopes
  if(inherits(isoscape, "isoStack")){
    ni = length(isoscape)
    
    for(i in 1:ni){
      if(inherits(isoscape[[i]], c("RasterStack", "RasterBrick"))){
        warning("raster objects are depreciated, transition to package terra")
        if(inherits(i, "RasterStack")){
          crs(isoscape[[i]][[1]]) = crs(isoscape[[i]][[2]]) = 
            crs(isoscape[[i]])
        }
        isoscape[[i]] = rast(isoscape[[i]])
      } 
      if(nlyr(isoscape[[i]]) != 2) {
        stop("Input isoscapes should be SpatRaster with two layers 
         (mean and standard deviation)")
      }
    }
  } else{
    ni = 1
    
    if(inherits(isoscape, "rescale")){
      isoscape = isoscape$isoscape.rescale
    }
    
    if(inherits(isoscape, c("RasterStack", "RasterBrick"))) {
      warning("raster objects are depreciated, transition to package terra")
      isoscape = rast(isoscape)
    } 
    
    if(inherits(isoscape, "SpatRaster")){
      if(is.na(crs(isoscape))) {
        stop("isoscape must have valid coordinate reference system")
      }
      if(nlyr(isoscape) != 2){
        stop("Input isoscape should be SpatRaster with two layers 
         (mean and standard deviation)")
      }
    } else {
      stop("isoscape should be a SpatRaster")
    }
  }
  
  #check known for multi-isotope
  if(ni > 1){
    #two options for ni>1, list of SODs, check and unpack each to spdf
    if(inherits(known, "list")){
      if(length(known) != ni){
        stop("length of known must equal length of isoStack")
      }
      #convert each SOD into SPDF
      for(i in 1:ni){
        known[[i]] = withCallingHandlers(
          message = addm,
          warning = addw,
          check_SOD(known[[i]], isoscape[[i]], bySite)
        ) 
      }
      #merge by sample - ?support partial data?
      k.spdf = known[[1]]
      kmlen = length(known[[1]])
      for(i in 2:ni){
        k.spdf = merge(k.spdf, known[[i]]@data, by = "Sample_ID", 
                       all.x = FALSE)
        if(bySite){
          if(!all(k.spdf$Site_ID.x == k.spdf$Site_ID.y)){
            stop("different Site_ID values for same samples")
          }
          k.spdf = k.spdf[, names(k.spdf) != "Site_ID.x"]
          names(k.spdf)[names(k.spdf) == "Site_ID.y"] = "Site_ID"
        }
        #move Sample_IDs to last column
        k.spdf@data = cbind(k.spdf@data[,-1], k.spdf@data[1])
        kmlen = max(kmlen, length(known[[i]]))
      }
      if(length(k.spdf) != kmlen){
        warning(paste("non-matching samples in known,", length(k.spdf),
                      "of", kmlen, "samples being used"))
      }
      known = k.spdf
    } 
    #if ni == 1
  } else{
    if(inherits(known, "subOrigData")){
      known = withCallingHandlers(
        message = addm,
        warning = addw,
        check_SOD(known, isoscape, bySite)
      )       
    }
  }
  #SOD or SOD list will now be converted to this
  if(inherits(known, "SpatialPointsDataFrame")){
    known = withCallingHandlers(
      message = addm,
      warning = addw,
      if(ni > 1){
        check_SPDF(known, isoscape[[1]], bySite, ni)
      } else{
        check_SPDF(known, isoscape, bySite, ni)
      }
    )
  } else{
    stop("invalid object provided for known")
  }

  #check recal
  if(!inherits(recal, "logical")){
    stop("recal must be logical")
  }
  if(!recal){
    valiTime = nrow(known)
    valiStation = 1
    bySite = FALSE
  }
  
  #check valiTime
  if(valiTime < 2){
    stop("valiTime must be an integer greater than 1")
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
  
  #check by
  if(!(as.integer(by) == by) || by < 1 || by > 25){
    stop("by must be an integer between 1 and 25")
  }
  
  #check name
  if(!is.null(name)){
    if(!inherits(name, "character")){
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
  if(recal){
    val_stations = sort(sample(ids, valiStation, replace = FALSE))
    for (i in seq_len(valiTime)[-1]){
      val_stations = rbind(val_stations, 
                           sort(sample(ids, valiStation, 
                                       replace = FALSE)))
    }
  } else{
    val_stations = matrix(ids, nrow = rowLength)
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
    
    if(recal){
      if(ni > 1){
        rescales = list()
        for(j in 1:ni){
          m_sub = m
          m_sub@data = m_sub@data[,(j * 2 - 1):(j * 2)]
          class(m_sub) = "QAData"
          rescales[[j]] = withCallingHandlers(
            message = addm,
            warning = addw,
            calRaster(m_sub, isoscape[[j]], mask, genplot = FALSE, 
                      verboseLM = FALSE)[[1]]
          )
        }
        rescale = isoStack(rescales)
      } else{
        class(m) = "QAData"
        rescale = withCallingHandlers(
          message = addm,
          warning = addw,
          calRaster(m, isoscape, mask, genplot = FALSE, 
                    verboseLM = FALSE)
        )
      } 
    } else{
      rescale = isoscape
    }

    pd = withCallingHandlers(
      message = addm,
      warning = addw,
      pdRaster(rescale, unknown = 
                 data.frame(row.names(v@data), 
                            v@data[,seq(1, ni*2-1, by=2)]), 
               prior = prior, genplot = FALSE)
    )
    
    v = vect(v)

    # pd value for each validation sample or site
    pd_temp = double(nlyr(pd))
    for(j in seq_along(pd_temp)){
      pd_temp[j] = extract(pd[[j]], v[j], factors = FALSE)[,2]
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
      rv_temp = double(nlyr(qtl))
      pre_temp = double(nlyr(qtl))
      for(k in seq_len(nlyr(qtl))){
        rv_temp[k] = extract(qtl[[k]], v[k], 
                             factors = FALSE)[, 2]
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
        if(any(!is.na(rv_temp))){
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
      rv_temp = double(nlyr(qtl))
      for(k in seq_along(rv_temp)){
        rv_temp[k] = extract(qtl[[k]], v[k], factors = FALSE)[, 2]
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
        if(any(!is.na(rv_temp))){
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

  #clean up warnings and messages
  wstack = unique(wstack)
  mstack = unique(mstack)
  trsh = lapply(wstack, warning, call. = FALSE)
  cat("\n")
  message(mstack)
  
  random_prob_density=1/length(cells(isoscape[[1]]))

  result = list(name, val_stations, pd_v, prption_byArea, 
                prption_byProb, precision, random_prob_density, by)
  names(result) = c("name", "val_stations", "pd_val", 
                    "prption_byArea", "prption_byProb", "precision",
                    "random_prob_density", "by")
  class(result) = "QA"
  return(result)
}

check_SOD = function(known, isoscape, bySite){
  
  col_site = NULL
  
  #check quality of SOD
  if(is.null(known$data) || is.null(known$marker)){
    stop("missing information in subOrigData known")
  }
  #find data columns
  col_m = match(known$marker, names(known$data))
  col_sd = match(paste0(known$marker, ".sd"), names(known$data))
  if(is.na(col_m) | is.na(col_sd)){
    stop("cannot match marker to data table in subOrigData known")
  }
  #find site column
  if(bySite){
    col_site = match("Site_ID", names(known$data))
    if(is.na(col_site)){
      stop("no Site_ID field in known; provide Site_IDs or use bySite = FALSE")
    }
  }
  #find samples column
  col_sample = match("Sample_ID", names(known$data))
  #pull out SPDF
  known = known$data
  #simplify data slot
  known@data = known@data[, c(col_m, col_sd, col_site, col_sample)]
  #check projection
  if(is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } 
  if(proj4string(known) != crs(isoscape, proj = TRUE)){
    known = spTransform(known, crs(isoscape))
    message("known was reprojected")
  } 
  
  return(known)
}

check_SPDF = function(known, isoscape, bySite, ni){
  
  col_site = NULL
  
  #check for enough columns
  if(ncol(known@data) < ni * 2){
    if(ncol(known@data) == 1 & ni == 1){
      warning("use of known with 1 data column is depreciated; known 
            should include the measured 
            isotope values and 1 sd uncertainty for each
            isotope system; assuming equal uncertainty for all samples")
      known@data[,2] = rep(0.0001, nrow(known@data))
    } else{
      stop("known must include the measured 
            isotope values and 1 sd uncertainty for each sample")
    }
  }
  
  #check all data columns are numeric w/ no missing values
  for(i in 1:(ni*2)){
    if(!is.numeric(known@data[,i])){
      stop("non-numeric data in sample value fields of known")
    }
    if(any(is.na(known@data[, i])) || 
       any(is.nan(known@data[, i])) || 
       any(is.null(known@data[, i]))){
      stop("Missing values detected in known sample value fields")
    }
  }
  
  #check that all SD values are greater than zero
  for(i in seq(2, ni*2, by = 2)){
    if(any(!(known@data[, i] > 0))){
      stop("negative or zero values found in known uncertainties")
    }
  }
  
  #check for Site_ID column
  if(bySite){
    col_site = match("Site_ID", names(known@data))
    if(is.na(col_site)){
      stop("no Site_ID field in known; provide Site_IDs or use bySite = FALSE")
    }
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
  if(proj4string(known) != crs(isoscape, proj = TRUE)){
    known = spTransform(known, crs(isoscape))
    message("known was reprojected")
  } 
  if(nrow(known) < 10){
    warning("there are fewer than 10 known samples")
  }
  if(nrow(known) < 3){
    stop("QA requires at least 3 known samples")
  }
  
  return(known)
}
