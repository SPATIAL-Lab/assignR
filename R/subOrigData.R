subOrigData = function(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
                        age_code = NULL, mask = NULL) {
  
  #load data in funtion environ
  data("knownOrig_samples", envir = environment())
  data("knownOrig_sites", envir = environment())
  data("knownOrig_sources", envir = environment())
  knownOrig_samples = knownOrig_samples
  knownOrig_sites = knownOrig_sites
  knownOrig_sources = knownOrig_sources
  
  result = knownOrig_samples

  if(!marker %in% colnames(knownOrig_samples)){
    stop("marker must be column name for isotope data field")
  }
  
  if(!is.null(taxon)){
    if(!all(taxon %in% unique(knownOrig_samples$Taxon))){
      warning("One or more taxa not present in database")
    }
    result = knownOrig_samples[knownOrig_samples$Taxon %in% taxon,]
  }
  
  if(!is.null(group)){
    if(!all(group %in% unique(knownOrig_samples$Group))){
      warning("One or more groups not present in database")
    }
    r = knownOrig_samples[knownOrig_samples$Group %in% group,]
    result = result[result$Sample_ID %in% r$Sample_ID,] 
  }
  
  if(!is.null(reference)){
    if(!all(reference %in% unique(knownOrig_sources$Dataset_ID))){
      warning("One or more references not present in database")
    }
    result = result[result$Dataset_ID %in% reference,] 
  }
  
  if(!is.null(age_code)){
    if(!all(age_code %in% unique(knownOrig_samples$Age_class))){
      warning("One or more age codes not present in database")
    }
    r = knownOrig_samples[knownOrig_samples$Age_class %in% age_code,]
    result = result[result$Sample_ID %in% r$Sample_ID,]
  }
  
  result = result[!is.na(result[,marker]),]
  if(length(result) == 0){
    stop("No samples match query")
  }

  if(!is.null(mask)) {
    if(class(mask)[1] == "SpatialPolygonsDataFrame" || 
       class(mask)[1] == "SpatialPolygons"){
      if(is.na(proj4string(mask))){
        stop("mask must have coordinate reference system")
      } else if(!identicalCRS(knownOrig_sites, mask)) {
        mask = spTransform(mask, proj4string(knownOrig_sites))
        warning("mask was reprojected")
      }
      result_sites = knownOrig_sites[mask,]
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }

    if(length(result_sites) > 0) {
      if(is.null(result)){
        result = knownOrig_samples[knownOrig_samples$Site_ID %in% 
                                     result_sites$Site_ID,]
      } else{
        result = result[result$Site_ID %in% result_sites$Site_ID,]
      }
      
      if(length(result) > 0) {
        plot(mask, axes = TRUE)
        plot(result_sites, add = TRUE, col = "red")
      } else {
        stop("No samples found in mask\n")
      }
    } else {
      stop("No sites found in mask\n")
    }
    
  } else {
    result_sites = knownOrig_sites[knownOrig_sites$Site_ID %in%
                                     result$Site_ID,]
    plot(wrld_simpl, axes = TRUE)
    plot(result_sites, add = TRUE, col = "red", cex = 0.5)
    
  }
  
  message(paste(length(result[,1]),"samples are found from", 
          length(result_sites), "sites"))
  
  result_data = merge(result_sites, result, by = "Site_ID", 
                      all.x = FALSE, duplicateGeoms = TRUE)
  result_sources = knownOrig_sources[knownOrig_sources$Dataset_ID %in%
                                       result_data$Dataset_ID,]
  
  return(list(result_data, result_sources))
}
