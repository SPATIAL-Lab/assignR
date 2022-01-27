wDist = function(pdR, sites){
  
  if(class(pdR) != "RasterLayer" & class(pdR) != "RasterStack" & class(pdR) != "RasterBrick"){
    stop("input probability density map (pdR) should be one of the following class: RasterLayer, RasterStack or RasterBrick")
  }
  if(is.na(proj4string(pdR))){
    stop("pdR must have coord. ref.")
  }
  
  if(class(sites)[1] == "SpatialPoints" || 
     class(sites)[1] == "SpatialPointsDataFrame"){
    if(length(sites) != nlayers(pdR)){
      stop("sites and pdR have different lenghts; wDist requires one site per pdR layer")
    }
    if(is.na(proj4string(sites))){
      stop("sites must have coord. ref.")
    }
    if(proj4string(sites) != proj4string(pdR)){
      sites = spTransform(sites, crs(pdR))
    }
  } else{
    stop("sites should be a SpatialPoints object")
  }
  
  #make space
  wd = list() 
  
  p = function(y, w){
    Position(function(x) x >= y, w)
  }
  
  for(i in seq_along(sites)){
    pdSP = rasterToPoints(pdR[[i]], spatial = TRUE)
    
    #for safety, using projected data works on most platforms
    pdSP = spTransform(pdSP, "+proj=longlat +ellps=WGS84")
    sites = spTransform(sites, "+proj=longlat +ellps=WGS84")
    
    d = distGeo(pdSP, sites[i,])
    b = bearing(pdSP, sites[i,])
    w = pdSP@data[,1]
    d.dens = density(d, weights = w)
    b.dens = density(b, weights = w)    
    
    #record weighted mean of distance distribution
    s = weighted.mean(d, w)
    
    #find and record quantiles within weighted distance distribution
    dw = cbind(d, w)
    dw = dw[order(d),]
    for(j in 2:nrow(dw)){
      dw[j, 2] = dw[j, 2] + dw[j-1, 2]
    }

    qts = sapply(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), p, w = dw[, 2])
    s = c(s, dw[qts, 1])
    
    bw = cbind(b, w)

    #find minimum weight value in bearing data to establish 'break'
    bbins = seq(-180, 170, by = 10)
    mp = 1
    for(j in bbins){
      pbin = sum(bw[bw[,1] >= j & bw[,1] < j + 10, 2])
      if(pbin < mp){
        mp = pbin
        mpb = j
      }
    }
    
    #re-reference bearing data to break
    bw[, 1] = bw[, 1] - mpb
    for(j in seq_along(bw[,1])){
      if(bw[j, 1] < 0){
        bw[j, 1] = bw[j, 1] + 360
      }
    }
    
    #weighted mean, re-referenced
    s = c(s, weighted.mean(bw[, 1], bw[, 2]))
    
    #find and record quantiles within weighted bearing distribution
    bw = bw[order(bw[, 1]),]
    for(j in 2:nrow(bw)){
      bw[j, 2] = bw[j, 2] + bw[j-1, 2]
    }
    qts = sapply(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), p, w = bw[, 2])
    s = c(s, bw[qts, 1])
    
    #rectify bearings
    s[9:16] = s[9:16] + mpb
    for(j in 9:16){
      if(s[j] >= 180){
        s[j] = s[j] - 360
      }
    }
    names(s) = c("wMeanDist", "w05Dist", "w10Dist", "w25Dist",
                  "w50Dist", "w75Dist", "w90Dist", "w95Dist", "wMeanBear",
                  "w05Bear", "w10Bear", "w25Bear", "w50Bear", "w75Bear",
                  "w90Bear", "w95Bear")
    
    wd[[i]] = list(stats = s, d.dens = d.dens, b.dens = b.dens)
  }
  
  class(wd) = "wDist"
  names(wd) = names(pdR)
  return(wd)
}

c.wDist = function(...){
  
  a = list(...)
  
  if(class(a[[1]])[1] != "wDist"){
    stop("... must be one or more wDist objects")
  }
  
  n = 0
  for(i in seq_len(length(a))){
    if(class(a[[i]])[1] == "wDist"){
      n = n + 1
    } else{
      stop("this method only accepts wDist objects as arguments")
    }
  }
  
  s = matrix(ncol = 16)
  k = 1
  for(i in seq_len(n)){
    nn = length(a[[i]])
    for(j in seq_len(length(a[[i]]))){
      if(k == 1){
        s = matrix(a[[i]][[j]]$stats, nrow = 1)
        nms = names(a[[i]][j])
      } else{
        s = rbind(s, a[[i]][[j]]$stats)
        nms = append(nms, names(a[[i]][j]))
      }
      k = k + 1
    }
  }
  
  s = as.data.frame(s)
  names(s) = names(a[[1]][[1]]$stats)
  s = cbind("Sample_ID" = nms, s)
  
  return(s)
}

plot.wDist = function(x, ..., bin = 20, pty = "both", index = c(1:5)){
  
  if(class(x) != "wDist"){
    stop("x must be a wDist object")
  }
  
  n = length(x)
  if(n == 0){
    stop("x is empty")
  }
  
  if(any(round(index) != index)){
    stop("index values must be integers")
  }
  if(length(index) > 5){
    message("more than 5 values in index, only the first 5 will be plotted")
    index = index[1:5]
  }
  if(length(index) == 5){
    if(all(index == c(1:5)) & n < 5){
      index = c(1:n)
    }
  }
  if(any(index > n)){
    message("index values exceeding length of x will not be plotted")
    index = index[index <= n]
  }
  
  if(!is.numeric(bin)){
    stop("bin must be numeric")
  }
  if(length(bin) > 1){
    stop("bin must be length 1")
  }
  if(bin <=0 | bin > 90){
    stop("bin must be a value between 0 and 90")
  }
  if(360 %% bin != 0){
    stop("bin should be a factor of 360")
  }
  
  if(!(pty %in% c("both", "dist", "bear"))){
    stop("pty not valid for plot.xist")
  }
  
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  
  if(pty %in% c("both", "dist")){
    #Distance
    d.xmax = d.ymax = 0
    d.dens = list()
    for(i in index){
      d.dens[[i]] = x[[i]]$d.dens
      d.xmax = max(d.xmax, max(d.dens[[i]]$x))
      d.ymax = max(d.ymax, max(d.dens[[i]]$y))
    }
    
    plot(d.dens[[index[1]]], xlim = c(0, d.xmax), ylim = c(0, d.ymax),
         main = "", ylab = "Probability density", xlab = "Distance (m)")
    for(i in index[-1]){
      lines(d.dens[[i]], col = i+1, )
    }
    legend("topright", legend = names(x), lty = 1, col = seq_len(n), 
           inset = 0.01)    
  }

  if(pty %in% c("both", "bear")){
    #Bearing
    b.dens = list()
    for(i in index){
      b.dens[[i]] = x[[i]]$b.dens
    }
    
    arc = function(a1, a2, b){
      a = seq(a1, a2, by = 0.5)
      r = 2 * pi * a / 360
      x = sin(r) * b
      y = cos(r) * b
      return(cbind(x, y))
    }
    
    wedge = function(a1, a2, b){
      xy = arc(a1, a2, b)
      xy = rbind(c(0,0), xy, c(0,0))
      return(xy)
    }
    
    bins = seq(-180, 179.9, by = bin)
    vals = numeric(length(bins))
    
    if(n > 3){
      mfr = 2
      if(n == 5){
        mfc = 3
      }else{
        mfc = 2
      }
    } else{
      mfr = 1
      mfc = n
    }
    par(mfrow = c(mfr, mfc), mar = c(1,1,3,1))
    
    for(i in index){
      b = b.dens[[i]]$x
      for(j in seq_along(b)){
        if(b[j] < -180){
          b[j] = b[j] + 360
        } else if(b[j] >= 180){
          b[j] = b[j] - 360
        }
      }
      y = b.dens[[i]]$y
      for(j in seq_along(bins)){
        vals[j] = sum(y[b >= bins[j] & b < bins[j] + bin])
      }
      
      b.max = max(vals)
      xy = arc(-180, 180, b.max)
      plot(xy, type = "l", col = "dark grey", axes = FALSE,
           xlab = "", ylab = "", asp = 1, main = names(x)[i])
      text(0.71 * b.max, 0.71 * b.max, signif(b.max, 2), 
           col = "dark grey", pos = 4, offset = 2)
      lines(arc(-180, 180, b.max/2), col = "dark grey")
      for(j in c(-180, -90, 0, 90)){
        lines(wedge(j, j, b.max * 1.05), col = "dark grey")
      }
      for(j in seq_along(bins)){
        xy = wedge(bins[j], bins[j] + bin, vals[j])
        c = col2rgb(i)
        polygon(xy, col = rgb(c[1], c[2], c[3],
                              alpha = 200, maxColorValue = 255))
      }
    }
  }

  return()
}
