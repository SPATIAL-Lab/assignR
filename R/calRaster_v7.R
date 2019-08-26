calRaster <- function (known, isoscape, mask = NULL, interpMethod = 2,
          NA.value = NA, ignore.NA = TRUE, genplot = TRUE, savePDF = FALSE, verboseLM = TRUE)
{
  #check that isoscape is valid and has defined CRS
  if(class(isoscape) == "RasterStack" | class(isoscape) == "RasterBrick") {
    if(is.na(proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
  } else {
    stop("isoscape should be a RasterStack or RastrBrick")
  }

  #check that known is valid and has defined, correct CRS
  if(class(known) != "SpatialPointsDataFrame") {
    stop("known should be a SpatialPointsDataFrame, see help page of calRaster function")
  }
  if(is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } 
  if(proj4string(known) != proj4string(isoscape)){
    known = spTransform(known, crs(isoscape))
    warning("known was reprojected")
  } 
  if(ncol(known@data) != 1){
    stop("known must include a 1-column data frame containing only the isotope values")
  }

  #check that mask is valid and has defined, correct CRS
  if(!is.null(mask)) {
    if(class(mask) == "SpatialPolygonsDataFrame" || class(mask) == "SpatialPolygons"){
      if(is.na(proj4string(mask))) {
        stop("mask must have valid coordinate reference system")
      }
      if(proj4string(mask) != proj4string(isoscape)){
        mask = spTransform(mask, crs(isoscape))
        warning("mask was reprojected")
      }
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }
  }

  #check that other inputs are valid
  if(!interpMethod %in% c(1,2)){
    stop("interpMethod should be 1 or 2")
  }
  if(class(genplot) != "logical") {
    stop("genplot should be logical (T or F)")
  }
  if (class(savePDF) != "logical") {
    stop("savePDF should be logical (T or F)")
  }

  #check and set isoscape NA value if necessary
  if(!is.na(NA.value)) {
    values(isoscape)[values(isoscape) == NA.value] <- NA
  }

  #get dimensions
  nSample <- nrow(known)

  #create space for regression variables
  tissue.iso <- vector("numeric", length = nSample)
  isoscape.iso <- vector("numeric", length = nSample)
  null.iso <- NULL

  #populate the dependent variable values
  tissue.iso <- known@data[, 1]

  #populate the independent variable values
  if (interpMethod == 1) {
    isoscape.iso <- raster::extract(isoscape, known,
                                    method = "simple")
  } else {
    isoscape.iso <- raster::extract(isoscape, known,
                                    method = "bilinear")
  }

  #warn if some known site have NA isoscape values
  if (any(is.na(isoscape.iso[, 1]))) {
    na <- which(is.na(isoscape.iso[, 1]))
    cat("\n\n----------------------------------------------------------------\n")
    cat("Warning: NO data are found at following locations:\n")
    print(known@coords[na])
    if (!ignore.NA) {
      stop("Delete these data in known origin data or use a different isoscape that has values at these locations")
    }

    #remove na values before continuing
    isoscape.iso <- isoscape.iso[-na, ]
    tissue.iso <- tissue.iso[-na, ]
  }

  #fit the regression model
  lmResult <- lm(tissue.iso ~ isoscape.iso[, 1])

  #output
  if (verboseLM){
    cat("\n\n---------------------------------------------------------------------------------\n")
    cat("rescale function uses linear regression model, the summary of this model is:\n")
    cat("---------------------------------------------------------------------------------\n")
    print(summary(lmResult))
  }

  #pull slope and intercept
  intercept <- as.numeric(coef(lmResult)[1])
  slope <- as.numeric(coef(lmResult)[2])

  #create rescaled prediction isoscape
  isoscape.rescale = isoscape[[1]] * slope + intercept

  #create data object for return
  x = isoscape.iso[, 1]
  y = tissue.iso
  xy = data.frame(x, y)

  if (genplot == TRUE || savePDF == TRUE) {
    #formatted lm equation for plotting
    equation <- function(mod) {
      lm_coef <- list(a = as.numeric(round(coef(mod)[1], digits = 2)),
                      b = as.numeric(round(coef(mod)[2], digits = 2)), r2 = round(summary(mod)$r.squared,
                                                                                  digits = 2))
      lm_eq <- substitute(italic(y) == a + b %.% italic(x) *
                            "," ~ ~italic(R)^2 ~ "=" ~ r2, lm_coef)
      as.expression(lm_eq)
    }
    #coordinates for placing equation legend in plot
    xl = max(x)
    yl = min(y) + 0.05 * diff(range(y))
  }
  
  if(genplot == TRUE){  
    #plot
    plot(x, y, pch = 21, bg="grey", xlab="Isoscape value", ylab="Tissue value", main="Rescale regression model")
    abline(lmResult)
    text(xl, yl, equation(lmResult), pos=2)
  }

  #combine uncertainties of isoscape and rescaling function
  sd <- (isoscape[[2]]^2 + (summary(lmResult)$sigma)^2)^0.5

  #stack the output rasters and apply names
  isoscape.rescale <- stack(isoscape.rescale, sd)
  names(isoscape.rescale) <- c("mean", "sd")

  #crop output if required
  if (!is.null(mask)) {
    isoscape.rescale <- crop(isoscape.rescale, mask)
  }

  #plot the output rasters
  if (genplot == TRUE) {
    print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE),
                 main = "Rescaled mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "Rescaled sd"))
  }

  #pdf output
  if (savePDF == TRUE) {
    dir.create("output")
    pdf("./output/rescale.result.pdf", width = 6, height = 4)
    
    #plot
    plot(x, y, pch = 21, bg="grey", xlab="Isoscape value", ylab="Tissue value", main="Rescale regression model")
    abline(lmResult)
    text(xl, yl, equation(lmResult), pos=2)
    
    print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE),
                 main = "Rescaled mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "Rescaled sd"))
    dev.off()
  }

  #set names for return data object
  names(xy) <- c("isoscape.iso", "tissue.iso")

  #package results
  result <- list(isoscape.rescale = isoscape.rescale, lm.data = xy,
                lm.model = lmResult)
  class(result) <- "rescale"
          
  #done
  return(result)
}
