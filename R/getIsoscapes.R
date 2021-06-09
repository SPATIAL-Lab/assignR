getIsoscapes = function(isoType = "GlobalPrecipGS", timeout = 1200){

  dpath.pre = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/"
  
  if(!is.numeric(timeout)){
    stop("timeout must be a number")
  }
  
  if(!(isoType %in% names(GIconfig))){
   stop("isoType invalid")
  }
  
  giconfig = GIconfig[[match(isoType, names(GIconfig))]]

  wd = getwd()
  setwd(tempdir())
  ot = getOption("timeout")
  options(timeout = timeout)
  
  dlf = function(fp, fn, ot, wd){
    dfs = tryCatch({
      download.file(fp, fn)
    },
    error = function(cond){
      options(timeout = ot)
      setwd(wd)
      stop(cond)
    })
    return(dfs)
  }
  
  pdlf = function(dfs, wd, ot, isoType){
    if(dfs != 0){
      setwd(wd)
      options(timeout = ot)
      stop("Non-zero download exit status")
    }
  }
  
  if(!file.exists(giconfig$dpath.post)){
    dfs = dlf(paste0(dpath.pre, giconfig$dpath.post), giconfig$dpath.post, ot, wd)
    pdlf(dfs, wd, ot, isoType)
    dlflag = TRUE
  } else{
    dlflag = FALSE
  }
  
  procRest = function(fn, lnames, onames){
    if((!all(lnames %in% list.files())) | dlflag){
      uz = unzip(fn)
    }
    rs = list()
    for(i in 1:length(lnames)){
      rs[[i]] = raster(lnames[i])
    }  
    names(rs) = onames
    return(rs)
  }
  
  rs = tryCatch({
    procRest(giconfig$dpath.post, giconfig$lnames, giconfig$onames)    
  },
  error = function(cond){
    stop(cond)
  },
  finally = {
    options(timeout = ot)
    setwd(wd)
  })
  
  switch(giconfig$eType,
         { #1
           if(length(rs) > 1){
             out = stack(rs)
           } else{
             out = rs
           }
         },
         { #2
           out = stack(rs)
         })
  
  message(paste0("Refer to ", tempdir(), "\\metadata.txt for 
  documentation and citation information"))

  return(out)
}

