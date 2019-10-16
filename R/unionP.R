unionP <- function(pdR){
  if(class(pdR) != "RasterStack"){
    stop("input probability density map (pdR) should be RasterLayer")
  }
  result <- (1 - pdR[[1]])
  for(i in 2:raster::nlayers(pdR)){
    result <- raster::overlay(result, pdR[[i]], fun = function(x,y){return(x*(1-y))})
  }
  raster::plot(1-result)
  graphics::title("Union Probability")
  return(1-result)
}
