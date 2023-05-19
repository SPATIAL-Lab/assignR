internal_files = function(){
  #add data to package environment
  assign("wrld_simpl", terra::vect(system.file("extdata/wrld_simpl.shp", package = "assignR")),
         pos = "package:assignR")
  assign("states", terra::vect(system.file("extdata/states.shp", package = "assignR")),
         pos = "package:assignR")
  assign("naMap", terra::vect(system.file("extdata/naMap.shp", package = "assignR")),
         pos = "package:assignR")
  assign("d2h_lrNA", terra::rast(system.file("extdata/d2h_lrNA.tif", package = "assignR")),
         pos = "package:assignR")
  assign("sr_MI", terra::rast(system.file("extdata/sr_MI.tif", package = "assignR")),
         pos = "package:assignR")
  
  assign("knownOrig", list(sites = terra::vect(system.file("extdata/knownOrig_sites.shp", package = "assignR")), 
                         samples = read.csv(system.file("extdata/knownOrig_samples.csv", package = "assignR")), 
                         sources = read.csv(system.file("extdata/knownOrig_sources.csv", package = "assignR"))),
         pos = "package:assignR")
}

.onAttach = function(libname, pkgname){
  #run attach
  internal_files()
  
  #cleanup removes terra pointers from package environment at end of session 
  clean = function(e){
    detach("package:assignR", unload = TRUE, character.only = TRUE)
  }
  e = as.environment("package:assignR")
  g = function(){
    reg.finalizer(e, clean, onexit = TRUE)
  }
  g()
}

