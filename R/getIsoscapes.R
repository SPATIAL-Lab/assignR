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
  on.exit({
    options(timeout = ot)
    setwd(wd)
  })
  
  dlf = function(fp, fn){
    dfs = tryCatch({
      download.file(fp, fn)
    },
    warning = function(cond){
      return(cond)
    },
    error = function(cond){
      return(cond)
    })
    return(dfs)
  }
  
  pdlf = function(dfs, wd, ot){
    setwd(wd)
    options(timeout = ot)
    message(paste("Download failed with status/message: \n", dfs))
  }
  
  if(!file.exists(giconfig$dpath.post)){
    dfs = dlf(paste0(dpath.pre, giconfig$dpath.post), giconfig$dpath.post)
    if(!is.numeric(dfs)){
      pdlf(dfs, wd, ot)
      return(NULL)
    }else if(dfs != 0){
      pdlf(dfs, wd, ot)
      return(NULL)
    }
  }
  
  procRest = function(fn, lnames, onames){
    if(file.exists("zRec.txt")){
      zRec = readLines("zRec.txt")
    } else{
      zRec = "none"
    }
    if((!all(lnames %in% list.files())) | (zRec != fn)){
      uz = unzip(fn)
      writeLines(fn, "zRec.txt")
    }
    rs = list()
    for(i in 1:length(lnames)){
      rs[[i]] = rast(lnames[i])
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
             out = rast(rs)
           } else{
             out = rs
           }
         },
         { #2
           out = rast(rs)
         })
  
  message(paste0("Refer to ", tempdir(), "\\metadata.txt for 
  documentation and citation information"))
  
  return(out)
}
