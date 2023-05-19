internal_files = function(){
  assign("wrld_simpl", terra::vect(system.file("extdata/wrld_simpl.shp", package = "assignR")),
         envir = as.environment("package:assignR"))
  assign("states", terra::vect(system.file("extdata/states.shp", package = "assignR")),
         envir = as.environment("package:assignR"))
  assign("naMap", terra::vect(system.file("extdata/naMap.shp", package = "assignR")),
         envir = as.environment("package:assignR"))
  assign("d2h_lrNA", terra::rast(system.file("extdata/d2h_lrNA.tif", package = "assignR")),
         envir = as.environment("package:assignR"))
  assign("sr_MI", terra::rast(system.file("extdata/sr_MI.tif", package = "assignR")),
         envir = as.environment("package:assignR"))
  
  assign("knownOrig", list(sites = terra::vect(system.file("extdata/knownOrig_sites.shp", package = "assignR")), 
                         samples = read.csv(system.file("extdata/knownOrig_samples.csv", package = "assignR")), 
                         sources = read.csv(system.file("extdata/knownOrig_sources.csv", package = "assignR"))),
         envir = as.environment("package:assignR"))
}

.onAttach = function(libname, pkgname){
  internal_files()
}
