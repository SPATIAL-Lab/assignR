jointP = function(pdR){
  if(class(pdR)[1] != "RasterStack" & class(pdR)[1] != "RasterBrick"){
    stop("input probability density map (pdR) should be RasterStack or 
         RasterBrick")
  }
  n = nlayers(pdR)
  result = pdR[[1]] * pdR[[2]]
  if(n > 2){
    for(i in seq_len(n)[-1:-2]){
      result = result * pdR[[i]]
    }
  }
  result = result / cellStats(result,sum)
  names(result) = "Joint_Probability"
  p = options("scipen")
  on.exit(options(scipen = p))
  options(scipen = -2)
  plot(result)
  title("Joint Probability")
  return(result)
}

unionP = function(pdR){
  if(class(pdR)[1] != "RasterStack" & class(pdR)[1] != "RasterBrick"){
    stop("input probability density map (pdR) should be RasterStack or 
         RasterBrick")
  }
  result = (1 - pdR[[1]])
  n = nlayers(pdR)
  for(i in seq_len(n)[-1]){
    result = overlay(result, pdR[[i]], fun = function(x,y){return(x*(1-y))})
  }
  plot(1-result)
  title("Union Probability")
  return(1-result)
}