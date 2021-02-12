calRaster = function (known, isoscape, mask = NULL, interpMethod = 2,
                      NA.value = NA, ignore.NA = TRUE, genplot = TRUE, 
                      outDir = NULL, verboseLM = TRUE){

  #check that isoscape is valid and has defined CRS
  if(class(isoscape)[1] == "RasterStack" | class(isoscape)[1] == "RasterBrick") {
    if(is.na(proj4string(isoscape))) {
      stop("isoscape must have valid coordinate reference system")
    }
    if(nlayers(isoscape) != 2) {
      stop("isoscape should be RasterStack or RasterBrick with two layers 
         (mean and standard deviation)")
    }
  } else {
    stop("isoscape should be a RasterStack or RasterBrick")
  }

  #check that known is valid and has defined, correct CRS
  if(!(class(known)[1] %in% c("subOrigData", "QAData", 
                              "SpatialPointsDataFrame"))) {
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
    known = known$data
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
    if(class(known) == "QAData"){
      class(known) = "SpatialPointsDataFrame"
    } else{
      warning("user-provided known; assuming measured isotope value and 1 sd
            uncertainty are contained in columns 1 and 2, respectively")
    }
    col_m = 1
    col_sd = 2
  }
  if(any(is.na(known@data[, col_m])) || 
     any(is.nan(known@data[, col_m])) || 
     any(is.null(known@data[, col_m]))){
    stop("Missing values detected in known values")
  }
  if(any(is.na(known@data[, col_sd])) || 
     any(is.nan(known@data[, col_sd])) || 
     any(is.null(known@data[, col_sd]))){
    stop("Missing values detected in known uncertainties")
  }
  if(any(known@data[, col_sd] == 0)){
    stop("zero values found in known uncertainties")
  }
  if(is.na(proj4string(known))) {
    stop("known must have valid coordinate reference system")
  } 
  if(proj4string(known) != proj4string(isoscape)){
    known = spTransform(known, crs(isoscape))
    warning("known was reprojected")
  } 

  #check that mask is valid and has defined, correct CRS
  if(!is.null(mask)) {
    if(class(mask)[1] == "SpatialPolygonsDataFrame" || 
       class(mask)[1] == "SpatialPolygons"){
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
  if(class(genplot)[1] != "logical") {
    stop("genplot should be logical (T or F)")
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
  
  #extract with mask
  if(!is.null(mask)){
    known = known[mask,]
    isoscape = crop(isoscape, mask)
  }

  #check and set isoscape NA value if necessary
  if(!is.na(NA.value)) {
    tempVals = getValues(isoscape)
    tempVals[tempVals == NA.value] = NA
    isoscape = setValues(isoscape, tempVals)
  }

  #get dimensions
  nSample = nrow(known)

  #create space for regression variables
  null.iso = NULL

  #populate the dependent variable values
  tissue.iso = known@data[, col_m]
  tissue.iso.wt = 1 / known@data[, col_sd]^2

  #populate the independent variable values
  if (interpMethod == 1) {
    isoscape.iso = extract(isoscape, known,
                                    method = "simple")
  } else {
    isoscape.iso = extract(isoscape, known,
                                    method = "bilinear")
  }

  #warn if some known sites have NA isoscape values
  if (any(is.na(isoscape.iso[, 1]))) {
    na = which(is.na(isoscape.iso[, 1]))
    wtxt = "NO isoscape values found at the following locations:\n"
    for(i in na){
      wtxt = paste0(wtxt, known@coords[i, 1], ", ", known@coords[i, 2], "\n")
    }
    if (ignore.NA) warning(wtxt)
    if (!ignore.NA) {
      cat(wtxt)
      stop("Delete these data in known origin data or use a different 
           isoscape that has values at these locations")
    }

    #remove na values before continuing
    tissue.iso = tissue.iso[!is.na(isoscape.iso[,1])]
    tissue.iso.wt = tissue.iso.wt[!is.na(isoscape.iso[,1])]
    isoscape.iso = isoscape.iso[!is.na(isoscape.iso[,1]), ]
    nSample = length(tissue.iso)
  }

  #fit the regression model
  lmResult = lm(tissue.iso ~ isoscape.iso[, 1], weights = tissue.iso.wt)

  #output
  if (verboseLM){
    cat("\n\n---------------------------------------
        ------------------------------------------\n")
    cat("rescale function uses linear regression model, 
        the summary of this model is:\n")
    cat("-------------------------------------------
        --------------------------------------\n")
    print(summary(lmResult))
  }

  #create data object for return
  x = isoscape.iso[, 1]
  y = tissue.iso
  w = tissue.iso.wt
  xyw = data.frame(x, y, w)

  if (genplot == TRUE || !is.null(outDir) ) {
    #formatted lm equation for plotting
    equation = function(mod) {
      lm_coef = list(a = as.numeric(round(coef(mod)[1], digits = 2)),
                      b = as.numeric(round(coef(mod)[2], digits = 2)), 
                      r2 = round(summary(mod)$r.squared, digits = 2))
      lm_eq = substitute(italic(y) == a + b %.% italic(x) *
                            "," ~ ~italic(R)^2 ~ "=" ~ r2, lm_coef)
      as.expression(lm_eq)
    }
    #coordinates for placing equation legend in plot
    xl = max(x)
    yl = min(y) + 0.05 * diff(range(y))
  }
  
  if(genplot == TRUE){  
    #plot
    plot(x, y, pch = 21, bg="grey", xlab="Isoscape value", 
         ylab="Tissue value", main="Rescale regression model")
    abline(lmResult)
    text(xl, yl, equation(lmResult), pos=2)
  }
  
  #pull slope and intercept
  intercept = as.numeric(coef(lmResult)[1])
  slope = as.numeric(coef(lmResult)[2])
  
  #create rescaled prediction isoscape
  isoscape.rescale = isoscape[[1]] * slope + intercept
  
  #simulate rescaling function variance
  isoscape.sim = matrix(0, nrow = nSample, ncol = 100)
  for(i in seq_along(isoscape.iso[,1])){
    isoscape.sim[i,] = rnorm(100, isoscape.iso[i, 1], isoscape.iso[i, 2])
  }
  isoscape.dev = tissue.dev = double()
  for(i in 1:100){
    lm.sim = lm(tissue.iso ~ isoscape.sim[,i], weights = tissue.iso.wt)
    isoscape.dev = c(isoscape.dev, isoscape.sim[,i] - isoscape.iso[,1])
    tissue.dev = c(tissue.dev, lm.sim$residuals)
  }

  #combine uncertainties of isoscape and rescaling function
  #rescaling variance is added variance from simulated fits
  sd = (isoscape[[2]]^2 + (var(tissue.dev) - var(isoscape.dev)))^0.5

  #stack the output rasters and apply names
  isoscape.rescale = stack(isoscape.rescale, sd)
  names(isoscape.rescale) = c("mean", "sd")

  #crop output if required
  if (!is.null(mask)) {
    isoscape.rescale = crop(isoscape.rescale, mask)
  }

  #plot the output rasters
  if (genplot == TRUE) {
    print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE),
                 main = "Rescaled mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "Rescaled sd"))
  }

  #pdf output
  if (!is.null(outDir)) {
    pdf(paste0(outDir, "/rescale_result.pdf"), width = 6, height = 4)
    
    #plot
    plot(x, y, pch = 21, bg="grey", xlab="Isoscape value", 
                   ylab="Tissue value", main="Rescale regression model")
    abline(lmResult)
    text(xl, yl, equation(lmResult), pos=2)
    
    print(spplot(isoscape.rescale$mean, scales = list(draw = TRUE),
                 main = "Rescaled mean"))
    print(spplot(isoscape.rescale$sd, scales = list(draw = TRUE),
                 main = "Rescaled sd"))
    dev.off()
  }

  #set names for return data object
  names(xyw) = c("isoscape.iso", "tissue.iso", "tissue.iso.wt")

  #package results
  result = list(isoscape.rescale = isoscape.rescale, lm.data = xyw,
                lm.model = lmResult)
  class(result) = c("rescale")
          
  #done
  return(result)
}
