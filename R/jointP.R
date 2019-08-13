jointP <- function(pdR){
  if(class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be RasterStack or RasterBrick")
  }
  n <- nlayers(pdR)
  result <- pdR[[1]]*pdR[[2]]
  if(n > 2){
    for(i in 3:n){
      result <- result*pdR[[i]]
    }
  }
  result <- result/cellStats(result,sum)
  names(result) <- "Joint_Probability"
  plot(result)
  title("Joint Probability")
  return(result)
}
