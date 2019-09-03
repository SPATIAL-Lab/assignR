pdRaster <- function(r, unknown, prior = NULL, mask = NULL, genplot = TRUE, saveFiles = FALSE) {
  if(class(r) == "rescale"){
    r <- r$isoscape.rescale
  }
  
  if(class(r) != "RasterStack" & class(r) != "RasterBrick"){
    stop("Input isoscape should be RasterStack or RasterBrick with two layers (mean and standard deviation)")
  } 
  if(nlayers(r) != 2) {
    stop("Input isoscape should be RasterStack or RasterBrick with two layers (mean and standard deviation)")
  }
  
  if(class(unknown) != "data.frame"){
    stop("unknown must be a data.frame")
  }
  if(any(is.na(unknown[,2])) || any(is.nan(unknown[,2])) || any(is.null(unknown[,2]))){
    stop("Missing isotope values detected in known")
  }
  
  if(!is.null(prior)){
    if(class(prior) != "RasterLayer"){
      stop("prior should be a raster with one layer")
    } 
    if(proj4string(prior) != proj4string(r[[1]])) {
      prior <- projectRaster(prior, crs=crs(r[[1]]))
      warning("prior was reprojected")
    }
    compareRaster(r[[1]], prior)
  }
  
  if(class(genplot) != "logical"){
    stop("genplot should be logical (T or F)")
  }
  
  if(class(saveFiles) != "logical"){
    stop("saveFiles should be logical (T or F)")
  }
  
  rescaled.mean = r[[1]]
  rescaled.sd = r[[2]]
  
  if (!is.null(mask)) {
    if(class(mask) == "SpatialPolygonsDataFrame" || class(mask) == "SpatialPolygons"){
      if (is.na(proj4string(mask))){
        stop("mask must have coord. ref.")
      } 
      if(proj4string(mask) != proj4string(r[[1]])) {
        mask <- spTransform(mask, crs(r[[1]]))
        warning("mask was reprojected")
      }
      rescaled.mean <- crop(rescaled.mean, mask)
      rescaled.sd <- crop(rescaled.sd, mask)
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
  
  if(saveFiles == TRUE){
    dir.create("output")
    dir.create("output/pdRaster_Gtif")
  }
  
  errorV <- getValues(rescaled.sd)
  meanV <- getValues(rescaled.mean)
  result <- NULL
  temp <- list()
  
  for (i in 1:n) {
    indv.data <- data[i, ]
    indv.id <- indv.data[1, 1]
    assign <- dnorm(indv.data[1, 2], mean = meanV, sd = errorV)
    if(!is.null(prior)){
      assign <- assign*getValues(prior)
    }
    assign.norm <- assign/sum(assign[!is.na(assign)])
    assign.norm <- setValues(rescaled.mean,assign.norm)
    if (i == 1){
      result <- assign.norm
    } else {
      result <- stack(result, assign.norm)
    }
    if(saveFiles == TRUE){
      filename <- paste("output/pdRaster_Gtif/", indv.id, ".like", ".tif", sep = "")
      writeRaster(assign.norm, filename = filename, format = "GTiff",
                  overwrite = TRUE)
    }
  }
  names(result) <- data[,1]

  if(saveFiles == TRUE){
    if (n > 5){
      pdf("./output/output_pdRaster.pdf", width = 10, height = 10)
      par(mfrow = c(ceiling(n/5), 5))
    } else {
      pdf("./output/output_pdRaster.pdf", width = 10, height = 10)
    }
  }

  if (genplot == TRUE){
    if (n == 1){
      pp <- spplot(result)
      print(pp)
    } else {
      for (i in 1:n){
        print(spplot(result@layers[[i]], scales = list(draw = TRUE), main=paste("Probability Density Surface for", data[i,1])))
      }
    }
  }
  
  if(saveFiles == TRUE){
    dev.off()
  }

  return(result)
}
