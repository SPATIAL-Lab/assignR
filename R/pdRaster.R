pdRaster <- function(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, outDir = NULL) {
  if(class(r) == "rescale"){
    r <- r$isoscape.rescale
  }
  
  if(class(r) != "RasterStack" & class(r) != "RasterBrick"){
    stop("Input isoscape should be RasterStack or RasterBrick with two layers (mean and standard deviation)")
  } 
  if(raster::nlayers(r) != 2) {
    stop("Input isoscape should be RasterStack or RasterBrick with two layers (mean and standard deviation)")
  }
  
  if(class(unknown) != "data.frame"){
    stop("unknown must be a data.frame")
  }
  if(any(is.na(unknown[,2])) || any(is.nan(unknown[,2])) || any(is.null(unknown[,2]))){
    stop("Missing isotope values detected in unknown")
  }
  
  if(!is.null(prior)){
    if(class(prior) != "RasterLayer"){
      stop("prior should be a raster with one layer")
    } 
    if(sp::proj4string(prior) != sp::proj4string(r[[1]])) {
      prior <- raster::projectRaster(prior, crs = raster::crs(r[[1]]))
      warning("prior was reprojected")
    }
    raster::compareRaster(r[[1]], prior)
  }
  
  if(class(genplot) != "logical"){
    stop("genplot should be logical (TRUE or FALSE)")
  }
  
  if(!is.null(outDir)){
    if(class(outDir) != "character"){
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
    if(class(mask) == "SpatialPolygonsDataFrame" || class(mask) == "SpatialPolygons"){
      if (is.na(sp::proj4string(mask))){
        stop("mask must have coord. ref.")
      } 
      if(sp::proj4string(mask) != sp::proj4string(r[[1]])) {
        mask <- sp::spTransform(mask, raster::crs(r[[1]]))
        warning("mask was reprojected")
      }
      rescaled.mean <- raster::crop(rescaled.mean, mask)
      rescaled.sd <- raster::crop(rescaled.sd, mask)
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }
  
  if (class(unknown) != "data.frame") {
    stop("unknown should be a data.frame, see help page of pdRaster function")
  }
  
  if (class(unknown) == "data.frame"){
    if(ncol(unknown) == 2){
      if(is.numeric(unknown[,2])){
        data <- unknown
      } else {
        stop("unknown column 2 must contain numeric values")
      }
    } else {
      stop("Wrong number of columns in unknown")
    }

  }
  
  n <- nrow(data)
  
  errorV <- raster::getValues(rescaled.sd)
  meanV <- raster::getValues(rescaled.mean)
  result <- NULL
  temp <- list()
  
  for (i in seq_len(n)) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, 1]
    assign <- stats::dnorm(indv.data[1, 2], mean = meanV, sd = errorV)
    if(!is.null(prior)){
      assign <- assign * raster::getValues(prior)
    }
    assign.norm <- assign / sum(assign[!is.na(assign)])
    assign.norm <- raster::setValues(rescaled.mean, assign.norm)
    if (i == 1){
      result <- assign.norm
    } else {
      result <- raster::stack(result, assign.norm)
    }
    if(!is.null(outDir)){
      filename <- paste0(outDir, "/", indv.id, "_like", ".tif", sep = "")
      raster::writeRaster(assign.norm, filename = filename, format = "GTiff", overwrite = TRUE)
    }
  }
  names(result) <- data[,1]

  if(!is.null(outDir)){
    if (n > 5){
      grDevices::pdf(paste0(outDir, "/output_pdRaster.pdf"), width = 10, height = 10)
      graphics::par(mfrow = c(ceiling(n/5), 5))
    } else {
      grDevices::pdf(paste0(outDir, "/output_pdRaster.pdf"), width = 10, height = 10)
    }
  }

  if (genplot == TRUE){
    if (n == 1){
      pp <- sp::spplot(result)
      print(pp)
    } else {
      for (i in seq_len(n)){
        print(sp::spplot(result@layers[[i]], scales = list(draw = TRUE), main=paste("Probability Density Surface for", data[i,1])))
      }
    }
  }
  
  if(!is.null(outDir)){
    grDevices::dev.off()
  }

  return(result)
}
