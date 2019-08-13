unionP <- function(pdR){
  if(class(pdR) != "RasterStack"){
    stop("input probability density map (pdR) should be RasterLayer")
  }
  result <- (1-pdR[[1]])
  for(i in 2:nlayers(pdR)){
    result <- overlay(result, pdR[[i]], fun = function(x,y){return(x*(1-y))})
  }
  plot(1-result)
  title("Union Probability")
  return(1-result)
}
