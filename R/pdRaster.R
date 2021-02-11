pdRaster = function(r, unknown, prior = NULL, mask = NULL, 
                    genplot = TRUE, outDir = NULL){
  UseMethod("pdRaster", r)
}

pdRaster.default = function(r, unknown, prior = NULL, mask = NULL, 
                           genplot = TRUE, outDir = NULL) {

  if(class(r)[1] == "rescale"){
    r = r$isoscape.rescale
  }
  
  if(class(r)[1] != "RasterStack" & class(r)[1] != "RasterBrick"){
    stop("Input isoscape should be RasterStack or RasterBrick with two layers 
         (mean and standard deviation)")
  } 
  if(nlayers(r) != 2) {
    stop("Input isoscape should be RasterStack or RasterBrick with two layers 
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
    rescaled.mean = crop(r[[1]], mask)
    rescaled.mean = mask(rescaled.mean, mask)
    rescaled.sd = crop(r[[2]], mask)
    rescaled.sd = mask(rescaled.sd, mask)
  }
  
  errorV = getValues(rescaled.sd)
  meanV = getValues(rescaled.mean)
  result = NULL
  temp = list()
  
  for (i in seq_len(n)) {
    indv.data = data[i, ]
    indv.id = indv.data[1,1]
    assign = dnorm(indv.data[1,2], mean = meanV, sd = errorV)
    if(!is.null(prior)){
      assign = assign * getValues(prior)
    }
    assign.norm = assign / sum(assign[!is.na(assign)])
    assign.norm = setValues(rescaled.mean, assign.norm)
    if (i == 1){
      result = assign.norm
    } else {
      result = stack(result, assign.norm)
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
  
  for(i in r){
    if(class(i)[1] != "RasterStack" & class(i)[1] != "RasterBrick"){
      stop("Input isoscapes should be RasterStack or RasterBrick with two layers 
         (mean and standard deviation)")
    } 
    if(nlayers(i) != 2) {
      stop("Input isoscapes should be RasterStack or RasterBrick with two layers 
         (mean and standard deviation)")
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
    rescaled.mean = crop(r[[1]][[1]], mask)
    rescaled.mean = mask(rescaled.mean, mask)
    rescaled.sd = crop(r[[1]][[2]], mask)
    rescaled.sd = mask(rescaled.sd, mask)
  }
  
  meanV = getValues(rescaled.mean)
  errorV = getValues(rescaled.sd)
  
  for(i in 2:ni){
    if(is.null(mask)){
      rescaled.mean = r[[i]][[1]]
      rescaled.sd = r[[i]][[2]]
    } else{
      rescaled.mean = crop(r[[i]][[1]], mask)
      rescaled.mean = mask(rescaled.mean, mask)
      rescaled.sd = crop(r[[i]][[2]], mask)
      rescaled.sd = mask(rescaled.sd, mask)
    }
    meanV = cbind(meanV, getValues(rescaled.mean))
    errorV = cbind(errorV, getValues(rescaled.sd))
  }

  result = NULL
  temp = list()
  assign = numeric(length(rescaled.mean))
  
  dev = cov(meanV, use = "pairwise.complete.obs")

  for (i in seq_len(n)) {
    indv.data = data[i, ]
    indv.id = indv.data[1, 1]
    indv.iso = indv.data[1, -1]
    
    for(j in seq_along(assign)){
      if(any(sapply(meanV[j, ], is.na))){
        assign[j] = NA
      } else{
        assign[j] = dmvn(as.numeric(indv.iso), meanV[j,], dev)
      }
    }

    if(!is.null(prior)){
      assign = assign * getValues(prior)
    }
    assign.norm = assign / sum(assign[!is.na(assign)])
    assign.norm = setValues(rescaled.mean, assign.norm)
    if (i == 1){
      result = assign.norm
    } else {
      result = stack(result, assign.norm)
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
  
  if(class(unknown)[1] == "refTrans"){
    unknown = process_refTrans(unknown)
  }
  
  if(class(unknown)[1] == "list"){
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
  
  if (class(unknown)[1] != "data.frame") {
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
    if(class(prior)[1] != "RasterLayer"){
      stop("prior should be a raster with one layer")
    } 
    if(proj4string(prior) != proj4string(r[[1]])) {
      prior = projectRaster(prior, crs = crs(r[[1]]))
      warning("prior was reprojected")
    }
    compareRaster(r[[1]], prior)
  }
  
  return(prior)
}

check_options = function(genplot, outDir){

  if(class(genplot)[1] != "logical"){
    stop("genplot should be logical (TRUE or FALSE)")
  }
  
  if(!is.null(outDir)){
    if(class(outDir)[1] != "character"){
      stop("outDir should be a character string")
    }
    if(!dir.exists(outDir)){
      warning("outDir does not exist, creating")
      dir.create(outDir)
    }
  }
  
  return()
}

check_mask = function(mask, r){

  if (!is.null(mask)) {
    if(class(mask)[1] == "SpatialPolygonsDataFrame" || 
       class(mask)[1] == "SpatialPolygons"){
      if (is.na(proj4string(mask))){
        stop("mask must have coord. ref.")
      } 
      if(proj4string(mask) != proj4string(r[[1]])) {
        mask = spTransform(mask, crs(r[[1]]))
        warning("mask was reprojected")
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
      pp = spplot(result)
      print(pp)
    } else {
      for (i in seq_len(n)){
        print(spplot(result@layers[[i]], scales = list(draw = TRUE), main=paste("Probability Density Surface for", data[i,1])))
      }
    }
  }
  
  if(!is.null(outDir)){
    dev.off()
  }
  
  return()
}

process_refTrans = function(unknown){
  if(class(unknown) == "refTrans"){
    if(ncol(unknown$data) == 3){
      un = data.frame("ID" = seq(1:nrow(unknown$data)), unknown$data[,1])
      names(un)[2] = names(unknown$data)[1]
      warning("no sample IDs in refTrans object, assigning numeric sequence")
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
