calRaster <- function (known, isoscape, mask = NULL, sdMethod = 2, interpMethod = 2,
          NA.value = NA, ignore.NA = TRUE, genplot = TRUE, savePDF = TRUE, verboseLM = T)
{
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

  #check that mask is valid and has defined, correct CRS
  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {
      if (is.na(proj4string(mask))) {
        stop("mask must have valid coordinate reference system")
      } else if(proj4string(mask) != proj4string(isoscape)){
        stop("mask must have same coordinate reference system as known")
      }
      known <- known[mask,]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }
  }

  #check that other inputs are valid
  if (sdMethod != 1 & sdMethod != 2) {
    stop("sdMethod should be 1 or 2. See help page for details")
  }
  if (class(genplot) != "logical") {
    stop("genplot should be logical (T or F)")
  }
  if (class(savePDF) != "logical") {
    stop("savePDF should be logical (T or F)")
  }

  #check and set isoscape NA value if necessary
  if (!is.na(NA.value)) {
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
  } else if (interpMethod == 2) {
    isoscape.iso <- raster::extract(isoscape, known,
                                    method = "bilinear")
  } else {
    stop("interpMethod should be either 1 or 2")
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

  if (genplot == TRUE) {
    #formatted lm equation for plotting
    equation <- function(mod) {
      lm_coef <- list(a = as.numeric(round(coef(mod)[1], digits = 2)),
                      b = as.numeric(round(coef(mod)[2], digits = 2)), r2 = round(summary(mod)$r.squared,
                                                                                  digits = 2))
      lm_eq <- substitute(italic(y) == a + b %.% italic(x) *
                            "," ~ ~italic(R)^2 ~ "=" ~ r2, lm_coef)
      as.character(as.expression(lm_eq))
    }

    #coordinate for placing equation legend in plot
    xmin = max(x) - 0.4 * (diff(range(x)))
    xmax = max(x)
    ymin = min(y) - 0.25 * diff(range(y))
    ymax = min(y) - 0.1 * diff(range(y))

    #plot object
    p11 <- ggplot(xy, aes(x = isoscape.iso[, 1], y = tissue.iso)) +
      geom_point(shape = 1) + geom_smooth(method = lm,
      se = T) + ggtitle("Rescale regression model") + theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(name = "Isoscape value") +
      scale_y_continuous(name = "Tissue value") + annotate("rect",
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = "white", colour = "red") + annotate("text",
      x = (xmax + xmin)/2, y = (ymax + ymin)/2, label = equation(lmResult),
      parse = TRUE)
    print(p11)
  }

  #interpret isoscape uncertainty values (sd or variance)
  if (sdMethod == 1) {
    sd1 <- isoscape[[2]]
  }
  if (sdMethod == 2) {
    sd1 <- isoscape[[2]]/2
  }

  #combine uncertainties of isoscape and rescaling function
  sd <- (sd1^2 + (summary(lmResult)$sigma)^2)^0.5

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
                 main = "rescale mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "rescale sd"))
  }

  #pdf output
  if (savePDF == TRUE) {
    dir.create("output")
    pdf("./output/rescale.result.pdf", width = 8, height = 4)
    print(p11)
    print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE),
                 main = "rescale mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "rescale sd"))
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
