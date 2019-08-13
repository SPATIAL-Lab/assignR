subOrigData <- function(taxon = NULL, group = NULL, mask = NULL) {
  data("knownOrig")
  result <- NULL
  if (!is.null(taxon)){
    if (all(taxon %in% unique(knownOrig$Taxon))) {
      result <- knownOrig[knownOrig$Taxon %in% taxon, ]
    } else {
      stop("taxon should be string or string vector given from Taxon column in knownOrig. Please see knownOrig help page!")
    }
  }
  if (!is.null(group)){
    if (all(group %in% unique(knownOrig$Group))) {
      result <- knownOrig[knownOrig$Group %in% group, ]
    } else {
      stop("group should be string or string vector given from Group column in knownOrig. Please see knownOrig help page!")
    }
  }

  if (!is.null(taxon) && !is.null(group))
    stop("Please either choose taxon or group")

  if (!is.null(mask)) {
    if (class(mask) == "SpatialPolygonsDataFrame") {
      if (is.na(proj4string(mask))){
        stop("mask must have coord. ref.")
      } else {
        mask <- spTransform(mask, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      }
      s <- data.frame(result@coords, result@data)
      o <- over(result, mask)
      overlap <- s[!is.na(o), ]
    } else {
      stop("mask should be a SpatialPolygonsDataFrame")
    }

    if (length(overlap[, 1]) != 0) {
      overlap <- SpatialPointsDataFrame(coords = cbind(overlap$Longitude,
                                                       overlap$Latitude), data = overlap, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      plot(mask, axes = T)
      plot(overlap, add = T, col = "red")
    } else {
      cat("No isotope data found in mask you choose!\n")
    }
    result <- overlap
  } else {
    require(maptools)
    data(wrld_simpl)
    plot(wrld_simpl, axes = T)
    points(result$Longitude, result$Latitude, col = "red", cex = 0.5)
  }
  cat(length(result[,1]),"data points are found\n")
  return(result[,3])
}
