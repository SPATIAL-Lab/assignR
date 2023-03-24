pdRaster = function(r, unknown, prior = NULL, mask = NULL, 
                    genplot = TRUE, outDir = NULL){
  UseMethod("pdRaster", r)
}

pdRaster.default = function(r, unknown, prior = NULL, mask = NULL, 
                           genplot = TRUE, outDir = NULL) {

  if(inherits(r, "rescale")){
    r = r$isoscape.rescale
  }
  
  ##legacy raster
  if(inherits(r, c("RasterStack", "RasterBrick"))) {
    warning("raster objects are depreciated, transition to package terra")
    r = rast(r)
    ##legacy raster
  } 
  
  if(inherits(r, "SpatRaster")){
    if(is.na(crs(r))){
      stop("r must have valid coordinate reference system")
    }
    if(nlyr(r) != 2) {
      stop("r should be a SpatRaster with two layers 
         (mean and standard deviation)")
    }
  } else{
    stop("r should be a SpatRaster with two layers 
         (mean and standard deviation)")
  }
  
  data = check_unknown(unknown, 1)
  n = nrow(data)
  
  prior = check_prior(prior, r)
  
  mask = check_mask(mask, r)  
  
  check_options(genplot, outDir)

  if(is.null(mask)){
    rescaled.mean = r[[1]]
    rescaled.sd = r[[2]]
  } else{
    rescaled.mean = crop(r[[1]], vect(mask))
    rescaled.mean = mask(rescaled.mean, vect(mask))
    rescaled.sd = crop(r[[2]], vect(mask))
    rescaled.sd = mask(rescaled.sd, vect(mask))
  }
  
  errorV = values(rescaled.sd, mat = FALSE)
  meanV = values(rescaled.mean, mat = FALSE)
  result = NULL
  temp = list()
  
  for (i in seq_len(n)) {
    indv.data = data[i, ]
    indv.id = indv.data[1,1]
    assign = dnorm(indv.data[1,2], mean = meanV, sd = errorV)
    if(!is.null(prior)){
      assign = assign * values(prior, mat = FALSE)
    }
    assign.norm = assign / sum(assign, na.rm = TRUE)
    assign.norm = setValues(rescaled.mean, assign.norm)
    if (i == 1){
      result = assign.norm
    } else {
      result = c(result, assign.norm)
    }
    if(!is.null(outDir)){
      filename = paste0(outDir, "/", indv.id, "_like", ".tif", sep = "")
      writeRaster(assign.norm, filename = filename, format = "GTiff", overwrite = TRUE)
    }
  }
  names(result) = data[,1]

  write_out(outDir, genplot, n, result, data)
  
  return(result)
}

pdRaster.isoStack = function(r, unknown, prior = NULL, mask = NULL, 
                            genplot = TRUE, outDir = NULL) {

  ni = length(r)
  
  for(i in seq(ni)){
    ##legacy raster
    if(inherits(r[[i]], c("RasterStack", "RasterBrick"))) {
      warning("raster objects are depreciated, transition to package terra")
      r[[i]] = rast(r[[i]])
      ##legacy raster
    }
    
    if(inherits(r[[i]], "SpatRaster")){
      if(is.na(crs(r[[i]]))){
        stop("isoscape must have valid coordinate reference system")
      }
      if(nlyr(r[[i]]) != 2) {
        stop("isoscape should be a SpatRaster with two layers 
         (mean and standard deviation)")
      }
    } else {
      stop("isoscape layers should be a SpatRaster")
    }
  }
  
  data = check_unknown(unknown, ni)
  n = nrow(data)
  
  prior = check_prior(prior, r[[1]])
  
  mask = check_mask(mask, r[[1]])  
  
  check_options(genplot, outDir)
  
  if(is.null(mask)){
    rescaled.mean = r[[1]][[1]]
    rescaled.sd = r[[1]][[2]]
  } else{
    rescaled.mean = crop(r[[1]][[1]], vect(mask))
    rescaled.mean = mask(rescaled.mean, vect(mask))
    rescaled.sd = crop(r[[1]][[2]], vect(mask))
    rescaled.sd = mask(rescaled.sd, vect(mask))
  }
  
  meanV = values(rescaled.mean, mat = FALSE)
  errorV = values(rescaled.sd, mat = FALSE)

  for(i in 2:ni){
    if(is.null(mask)){
      rescaled.mean = r[[i]][[1]]
      rescaled.sd = r[[i]][[2]]
    } else{
      rescaled.mean = crop(r[[i]][[1]], vect(mask))
      rescaled.mean = mask(rescaled.mean, vect(mask))
      rescaled.sd = crop(r[[i]][[2]], vect(mask))
      rescaled.sd = mask(rescaled.sd, vect(mask))
    }
    meanV = cbind(meanV, values(rescaled.mean, mat = FALSE))
    errorV = cbind(errorV, values(rescaled.sd, mat = FALSE))
  }

  result = NULL
  temp = list()
  assign = as.numeric(rep(NA, nrow(meanV)))
  cells = seq_along(meanV[,1])
  cellmask = apply(cbind(meanV, errorV), 1, anyNA)
  cells = cells[!cellmask]
  
  #sanity check
  cd = cdt = cor(meanV, use = "pairwise.complete.obs")^2
  diag(cdt) = NA
  if(any(cdt > 0.7, na.rm = TRUE)){
    warning("two or more isoscapes have shared variance > 0.7, added information
            will be limited, and specificity of assignments may be inflated")
  }
  
  dev = d.cell = cov(meanV, use = "pairwise.complete.obs")
  v = sqrt(diag(dev))
  d.l = list()
  for(i in cells){
    v.cell = errorV[i,] / v
    for(j in 1:ni){
      for(k in 1:ni){
        d.cell[j, k] = dev[j, k] * v.cell[j] * v.cell[k]
      }
    }
    d.l[[i]] = d.cell
  }
  
  for (i in seq_len(n)) {
    indv.data = data[i, ]
    indv.id = indv.data[1, 1]
    indv.iso = indv.data[1, -1]
    
    for(j in cells){
      assign[j] = dmvn(as.numeric(indv.iso), meanV[j,], d.l[[j]])
    }

    if(!is.null(prior)){
      assign = assign * values(prior, mat = FALSE)
    }
    assign.norm = assign / sum(assign[!is.na(assign)])
    assign.norm = setValues(rescaled.mean, assign.norm)
    if (i == 1){
      result = assign.norm
    } else {
      result = c(result, assign.norm)
    }
    if(!is.null(outDir)){
      filename = paste0(outDir, "/", indv.id, "_like", ".tif", sep = "")
      writeRaster(assign.norm, filename = filename, format = "GTiff", overwrite = TRUE)
    }
  }
  names(result) = data[,1]
  
  write_out(outDir, genplot, n, result, data)
  
  return(result)
}

check_unknown = function(unknown, n){
  
  if(inherits(unknown, "refTrans")){
    unknown = process_refTrans(unknown)
  }
  
  if(inherits(unknown, "list")){
    if(length(unknown) < n){
      stop("number of refTrans objects provided is less than number of isoscapes")
    }
    if(length(unknown) > n){
      warning(paste("more refTrans objects than isoscapes, only using first",
                    n, "objects"))
    }
    
    un = process_refTrans(unknown[[1]])
    nobs = nrow(un)
    
    for(i in 2:n){
      u = process_refTrans(unknown[[i]])
      if(nrow(u) != nobs){
        stop("different numbers of samples in refTrans objects")
      }
      un = merge(un, u, by.x = 1, by.y = 1)
      if(nrow(un) != nobs){
        stop("sample IDs in refTrans objects don't match")
      }
    }
    
    unknown = un
  }
  
  for(i in seq(ncol(unknown))){
    if(any(is.na(unknown[,i])) || any(is.nan(unknown[,i])) || 
       any(is.null(unknown[,i]))){
      stop("Missing values detected in unknown")
    }
  }
  
  if (!inherits(unknown, "data.frame")) {
    stop("unknown should be a data.frame, see help page of pdRaster function")
  }
  
  if(ncol(unknown) < (n + 1)){
    stop("unknown must contain sample ID in col 1 and marker values in col 2+")
  }
  if(ncol(unknown) > (n + 1)){
    warning("more than n marker + 1 cols in unknown, assuming IDs in col 1 and marker values in cols 2 to (n+1)")
  }
  for(i in 2:(n+1)){
    if(!is.numeric(unknown[,i])){
      stop("unknown data column(s) must contain numeric values")
    }
  }

  return(unknown[,1:(n+1)])
}

check_prior = function(prior, r){
  
  if(!is.null(prior)){
    #legacy raster
    if(inherits(prior, "RasterLayer")){
      warning("raster objects are depreciated, transition to package terra")
      prior = rast(prior)
      #legacy raster
    } 
    
    if(inherits(prior, "SpatRaster")){
      if(is.na(crs(prior))){
        stop("isoscape must have valid coordinate reference system")
      }
      if(crs(prior) != crs(r[[1]])) {
        prior = project(prior, crs = crs(r[[1]]))
        message("prior was reprojected")
      }
      compareGeom(prior, r)
    } else{
      stop("prior should be a SpatRaster")
    }
  }
  
  return(prior)
}

check_options = function(genplot, outDir){

  if(!inherits(genplot, "logical")){
    stop("genplot should be logical (TRUE or FALSE)")
  }
  
  if(!is.null(outDir)){
    if(class(outDir)[1] != "character"){
      stop("outDir should be a character string")
    }
    if(!dir.exists(outDir)){
      message("outDir does not exist, creating")
      dir.create(outDir)
    }
  }
  
  return()
}

check_mask = function(mask, r){

  if (!is.null(mask)) {
    if(inherits(mask, "SpatialPolygons")){
      if (is.na(proj4string(mask))){
        stop("mask must have valid coordinate reference system")
      } 
      if(proj4string(mask) != crs(r, proj = TRUE)){
        mask = spTransform(mask, crs(r))
        message("mask was reprojected")
      }
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }
  
  return(mask)
}

write_out = function(outDir, genplot, n, result, data){
  
  if(!is.null(outDir)){
    if (n > 5){
      pdf(paste0(outDir, "/output_pdRaster.pdf"), width = 10, height = 10)
      par(mfrow = c(ceiling(n/5), 5))
    } else {
      pdf(paste0(outDir, "/output_pdRaster.pdf"), width = 10, height = 10)
    }
  }
  
  if (genplot == TRUE){
    if (n == 1){
      pp = plot(result)
      print(pp)
    } else {
      pp = plot(result, main=paste("Probability Density Surface for", 
                                   data[, 1]))
      print(pp)
    }
  }
  
  if(!is.null(outDir)){
    dev.off()
  }
  
  return()
}

process_refTrans = function(unknown){
  if(inherits(unknown, "refTrans")){
    if(ncol(unknown$data) == 3){
      un = data.frame("ID" = seq(1:nrow(unknown$data)), unknown$data[,1])
      names(un)[2] = names(unknown$data)[1]
      message("no sample IDs in refTrans object, assigning numeric sequence")
    } else if(ncol(unknown$data) == 4){
      un = unknown$data[,1:2]
    } else{
      stop("incorrectly formatted refTrans object in unknown")
    }
  } else{
    stop("only refTrans objects can be provided in listed unknown")
  }
  return(un)
}
