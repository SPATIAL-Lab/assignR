calRaster = function (known, isoscape, mask = NULL, interpMethod = 2,
                      NA.value = NA, ignore.NA = TRUE, genplot = TRUE, 
                      outDir = NULL, verboseLM = TRUE){

  
  #check that isoscape is valid and has defined CRS
  ##legacy raster
  if(inherits(isoscape, c("RasterStack", "RasterBrick"))) {
    warning("raster objects are depreciated, transition to package terra")
    isoscape = rast(isoscape)
    ##legacy raster
  } 
  
  if(inherits(isoscape, "SpatRaster")){
    if(crs(isoscape) == ""){
      stop("isoscape must have valid coordinate reference system")
    }
    if(nlyr(isoscape) != 2) {
      stop("isoscape should be a SpatRaster with two layers 
         (mean and standard deviation)")
    }
  } else {
    stop("isoscape should be a SpatRaster")
  }

  #check that known is valid and has defined, correct CRS
  if(!(inherits(known, c("subOrigData", "QAData", 
                         "SpatialPointsDataFrame",
                         "SpatVector")))) {
    stop("known must be a subOrigData or SpatVector object")
  }
  if(inherits(known, "SpatialPointsDataFrame")){
    known = vect(known)
  }
  if(inherits(known, "subOrigData")){
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
    if(inherits(known, "QAData")){
      class(known) = "SpatVector"
    } else{
      message("user-provided known; assuming measured isotope value and 1 sd
            uncertainty are contained in columns 1 and 2, respectively")
    }
    if(ncol(known) < 2){
      if(ncol(known) == 1){
        warning("use of known with 1 data column is depreciated; known 
              should include a data frame containing the measured 
              isotope values (col 1) and 1 sd uncertainty (col 2); 
              assuming equal uncertainty for all samples")
        known[,2] = rep(0.0001, nrow(known))
      } else{
        stop("known must include a data frame containing the measured 
              isotope values (col 1) and 1 sd uncertainty (col 2)")
      }
    }
    if(any(!is.numeric(values(known)[,1]), !is.numeric(values(known)[,2]))){
      stop("known must include data containing the measured 
           isotope values (col 1) and 1 sd uncertainty (col 2)")
    }
    col_m = 1
    col_sd = 2
  }
  if(any(is.na(values(known)[, col_m])) || 
     any(is.nan(values(known)[, col_m])) || 
     any(is.null(values(known)[, col_m]))){
    stop("Missing values detected in known values")
  }
  if(any(is.na(values(known)[, col_sd])) || 
     any(is.nan(values(known)[, col_sd])) || 
     any(is.null(values(known)[, col_sd]))){
    stop("Missing values detected in known uncertainties")
  }
  if(any(values(known)[, col_sd] == 0)){
    stop("zero values found in known uncertainties")
  }
  if(crs(known) == "") {
    stop("known must have valid coordinate reference system")
  } 
  if(!same.crs(known, isoscape)){
    known = project(known, crs(isoscape))
    message("known was reprojected")
  } 

  #check that mask is valid and has defined, correct CRS
  mask = check_mask(mask, isoscape)

  #check that other inputs are valid
  if(!interpMethod %in% c(1,2)){
    stop("interpMethod should be 1 or 2")
  }
  if(!inherits(genplot, "logical")) {
    message("genplot should be logical (T or F), using default = T")
    genplot = TRUE
  }
  if(!is.null(outDir)){
    if(!inherits(outDir, "character")){
      stop("outDir should be a character string")
    }
    if(!dir.exists(outDir)){
      message("outDir does not exist, creating")
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
    tempVals = values(isoscape)
    tempVals[tempVals == NA.value] = NA
    isoscape = setValues(isoscape, tempVals)
  }

  #get dimensions
  nSample = nrow(known)

  #create space for regression variables
  null.iso = NULL

  #populate the dependent variable values
  tissue.iso = values(known)[, col_m]
  tissue.iso.sd = values(known)[, col_sd]
  tissue.iso.wt = 1 / tissue.iso.sd^2

  #populate the independent variable values
  if (interpMethod == 1) {
    isoscape.iso = extract(isoscape, known, method = "simple")[,2:3]
  } else {
    isoscape.iso = extract(isoscape, known, method = "bilinear")[,2:3]
  }
  #protect against negative values from interpolation
  isoscape.iso[,2] = pmax(isoscape.iso[,2], 
                          global(isoscape, min, na.rm = TRUE)[2, 1])

  #warn if some known sites have NA isoscape values
  if (any(is.na(isoscape.iso[, 1]))) {
    na = which(is.na(isoscape.iso[, 1]))
    wtxt = "NO isoscape values found at the following locations:\n"
    for(i in na){
      wtxt = paste0(wtxt, geom(known)[i, 3], ", ", geom(known)[i, 4], "\n")
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
  
  ti.corr = cor(isoscape.dev, tissue.dev)^2
  
  #combine uncertainties of isoscape and rescaling function
  #rescaling variance is frac of model variance uncorrelated w/ isoscape error
  rescale.sd = sqrt(isoscape[[2]]^2 + var(lmResult$residuals) * (1-ti.corr))
  
  #stack the output rasters and apply names
  isoscape.rescale = c(isoscape.rescale, rescale.sd)
  names(isoscape.rescale) = c("mean", "sd")

  #crop output if required
  if (!is.null(mask)) {
    isoscape.rescale = crop(isoscape.rescale, mask)
  }

  #plot the output rasters
  if (genplot == TRUE) {
    print(plot(isoscape.rescale, main = c("Rescaled mean", "Rescaled sd")))
  }

  #pdf output
  if (!is.null(outDir)) {
    pdf(paste0(outDir, "/rescale_result.pdf"), width = 6, height = 4)
    
    #plot
    plot(x, y, pch = 21, bg="grey", xlab="Isoscape value", 
                   ylab="Tissue value", main="Rescale regression model")
    abline(lmResult)
    text(xl, yl, equation(lmResult), pos=2)
    
    print(plot(isoscape.rescale, main = c("Rescaled mean", "Rescaled sd")))
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
