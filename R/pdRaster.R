pdRaster = function(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, outDir = NULL) {
  if(class(r) == "rescale"){
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
  
  if(any(is.na(unknown[,2])) || any(is.nan(unknown[,2])) || 
     any(is.null(unknown[,2]))){
    stop("Missing isotope values detected in unknown")
  }
  
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
  
  rescaled.mean = r[[1]]
  rescaled.sd = r[[2]]
  
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
      rescaled.mean = crop(rescaled.mean, mask)
      rescaled.sd = crop(rescaled.sd, mask)
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }
  
  if (class(unknown)[1] != "data.frame") {
    stop("unknown should be a data.frame, see help page of pdRaster function")
  }
  
  if (class(unknown)[1] == "data.frame"){
    if(ncol(unknown) == 2){
      if(is.numeric(unknown[,2])){
        data = unknown
      } else {
        stop("unknown column 2 must contain numeric values")
      }
    } else {
      stop("Wrong number of columns in unknown")
    }

  }
  
  n = nrow(data)
  
  errorV = getValues(rescaled.sd)
  meanV = getValues(rescaled.mean)
  result = NULL
  temp = list()
  
  for (i in seq_len(n)) {
    indv.data = data[i, ]
    indv.id = indv.data[1, 1]
    assign = dnorm(indv.data[1, 2], mean = meanV, sd = errorV)
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

  return(result)
}
