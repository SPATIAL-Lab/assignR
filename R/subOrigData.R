subOrigData <- function(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
                        age_code = NULL, mask = NULL) {
  
  result <- NULL
  
  if(!marker %in% colnames(knownOrig@data)){
    stop("marker must be column name for isotope data field")
  }
  
  if(!is.null(taxon)){
    if(!all(taxon %in% unique(knownOrig$Taxon))){
      warning("One or more taxa not present in database")
    }
    result = knownOrig[knownOrig$Taxon %in% taxon,]
  }
  if(!is.null(group)){
    if(!all(group %in% unique(knownOrig$Group))){
      warning("One or more groups not present in database")
    }
    r = knownOrig[knownOrig$Group %in% group,]
    if(is.null(result)){
      result = r
    }else{
      result = result[result$ID %in% r$ID,] 
    }
  }
  if(!is.null(reference)){
    if(!all(reference %in% unique(knownOrig$Reference))){
      warning("One or more references not present in database")
    }
    r = knownOrig[knownOrig$Reference %in% reference,]
    if(is.null(result)){
      result = r
    }else{
      result = result[result$ID %in% r$ID,] 
    }
  }
  if(!is.null(age_code)){
    if(!all(age_code %in% unique(knownOrig$Age_code))){
      warning("One or more age codes not present in database")
    }
    r = knownOrig[knownOrig$Age_code %in% age_code,]
    if(is.null(result)){
      result = r
    }else{
      result = result[result$ID %in% r$ID,] 
    }
  }
  
  result = result[!is.na(result@data[,marker]),]
  if(length(result) == 0){
    stop("No samples match query")
  }

  if(!is.null(mask)) {
    if(class(mask) == "SpatialPolygonsDataFrame" || class(mask) == "SpatialPolygons"){
      if(is.na(sp::proj4string(mask))){
        stop("mask must have coordinate reference system")
      } else {
        mask <- sp::spTransform(mask, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      }
      overlap <- result[mask,]
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }

    if(length(overlap) > 0) {
      sp::plot(mask, axes = TRUE)
      sp::plot(overlap, add = TRUE, col = "red")
    } else {
      stop("No samples found in mask\n")
    }
    
    result <- overlap
    
  } else {
    utils::data(wrld_simpl)
    sp::plot(wrld_simpl, axes = TRUE)
    sp::plot(result, add = TRUE, col = "red", cex = 0.5)
    
  }
  cat(length(result[,1]),"data points are found\n")
  return(result[,marker])
}
