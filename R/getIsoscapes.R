getIsoscapes = function(isoType = "GlobalPrecipGS", timeout = 1200){

  dpath.pre = "https://wateriso.utah.edu/waterisotopes/media/ArcGrids/"
  
  switch(isoType, 
         "GlobalPrecipGS" = {
           dpath.post = "GlobalPrecipGS.zip"
           lnames = c("d2h_GS.tif", "d2h_se_GS.tif", 
                      "d18o_GS.tif", "d18o_se_GS.tif")
           onames = c("d2h", "d2h.se", "d18o", "d18o.se")
           eType = 2
         },
         "USPrecipMA" =,
         "GlobalPrecipMA" = {
           dpath.post = switch(isoType,
                  "USPrecipMA" = "USPrecip.zip",
                  "GlobalPrecipMA" = "GlobalPrecip.zip")
           lnames = c("d2h_MA.tif", "d2h_se_MA.tif", 
                      "d18o_MA.tif", "d18o_se_MA.tif")
           onames = c("d2h_MA", "d2h_se_MA", "d18o_MA", "d18o_se_MA")
           eType = 2
         },
         "USPrecipMO" =,
         "GlobalPrecipMO" = {
           dpath.post = switch(isoType,
                               "USPrecipMO" = "USPrecip.zip",
                               "GlobalPrecipMO" = "GlobalPrecip.zip")
           lnames = c("d2h_01.tif", "d2h_se_01.tif", 
                      "d2h_02.tif", "d2h_se_02.tif",
                      "d2h_03.tif", "d2h_se_03.tif",
                      "d2h_04.tif", "d2h_se_04.tif",
                      "d2h_05.tif", "d2h_se_05.tif",
                      "d2h_06.tif", "d2h_se_06.tif",
                      "d2h_07.tif", "d2h_se_07.tif",
                      "d2h_08.tif", "d2h_se_08.tif",
                      "d2h_09.tif", "d2h_se_09.tif",
                      "d2h_10.tif", "d2h_se_10.tif",
                      "d2h_11.tif", "d2h_se_11.tif",
                      "d2h_12.tif", "d2h_se_12.tif",
                      "d18o_01.tif", "d18o_se_01.tif",
                      "d18o_02.tif", "d18o_se_02.tif",
                      "d18o_03.tif", "d18o_se_03.tif",
                      "d18o_04.tif", "d18o_se_04.tif",
                      "d18o_05.tif", "d18o_se_05.tif",
                      "d18o_06.tif", "d18o_se_06.tif",
                      "d18o_07.tif", "d18o_se_07.tif",
                      "d18o_08.tif", "d18o_se_08.tif",
                      "d18o_09.tif", "d18o_se_09.tif",
                      "d18o_10.tif", "d18o_se_10.tif",
                      "d18o_11.tif", "d18o_se_11.tif",
                      "d18o_12.tif", "d18o_se_12.tif")
           onames = c("d2h_01", "d2h_se_01", "d2h_02", "d2h_se_02", 
                      "d2h_03", "d2h_se_03", "d2h_04", "d2h_se_04",
                      "d2h_05", "d2h_se_05", "d2h_06", "d2h_se_06",
                      "d2h_07", "d2h_se_07", "d2h_08", "d2h_se_08",
                      "d2h_09", "d2h_se_09", "d2h_10", "d2h_se_10",
                      "d2h_11", "d2h_se_11", "d2h_12", "d2h_se_12",
                      "d18o_01", "d18o_se_01", "d18o_02", "d18o_se_02", 
                      "d18o_03", "d18o_se_03", "d18o_04", "d18o_se_04",
                      "d18o_05", "d18o_se_05", "d18o_06", "d18o_se_06",
                      "d18o_07", "d18o_se_07", "d18o_08", "d18o_se_08",
                      "d18o_09", "d18o_se_09", "d18o_10", "d18o_se_10",
                      "d18o_11", "d18o_se_11", "d18o_12", "d18o_se_12")
           eType = 2
         },
         "USPrecipALL" =,
         "GlobalPrecipALL" = {
           dpath.post = switch(isoType,
                               "USPrecipALL" = "USPrecip.zip",
                               "GlobalPrecipALL" = "GlobalPrecip.zip")
           lnames = c("d2h_MA.tif", "d2h_se_MA.tif",
                      "d2h_01.tif", "d2h_se_01.tif", 
                      "d2h_02.tif", "d2h_se_02.tif",
                      "d2h_03.tif", "d2h_se_03.tif",
                      "d2h_04.tif", "d2h_se_04.tif",
                      "d2h_05.tif", "d2h_se_05.tif",
                      "d2h_06.tif", "d2h_se_06.tif",
                      "d2h_07.tif", "d2h_se_07.tif",
                      "d2h_08.tif", "d2h_se_08.tif",
                      "d2h_09.tif", "d2h_se_09.tif",
                      "d2h_10.tif", "d2h_se_10.tif",
                      "d2h_11.tif", "d2h_se_11.tif",
                      "d2h_12.tif", "d2h_se_12.tif",
                      "d18o_MA.tif", "d18o_se_MA.tif",
                      "d18o_01.tif", "d18o_se_01.tif",
                      "d18o_02.tif", "d18o_se_02.tif",
                      "d18o_03.tif", "d18o_se_03.tif",
                      "d18o_04.tif", "d18o_se_04.tif",
                      "d18o_05.tif", "d18o_se_05.tif",
                      "d18o_06.tif", "d18o_se_06.tif",
                      "d18o_07.tif", "d18o_se_07.tif",
                      "d18o_08.tif", "d18o_se_08.tif",
                      "d18o_09.tif", "d18o_se_09.tif",
                      "d18o_10.tif", "d18o_se_10.tif",
                      "d18o_11.tif", "d18o_se_11.tif",
                      "d18o_12.tif", "d18o_se_12.tif")
           onames = c("d2h_MA", "d2h_se_MA", 
                      "d2h_01", "d2h_se_01", "d2h_02", "d2h_se_02", 
                      "d2h_03", "d2h_se_03", "d2h_04", "d2h_se_04",
                      "d2h_05", "d2h_se_05", "d2h_06", "d2h_se_06",
                      "d2h_07", "d2h_se_07", "d2h_08", "d2h_se_08",
                      "d2h_09", "d2h_se_09", "d2h_10", "d2h_se_10",
                      "d2h_11", "d2h_se_11", "d2h_12", "d2h_se_12",
                      "d18o_MA", "d18o_se_MA",
                      "d18o_01", "d18o_se_01", "d18o_02", "d18o_se_02", 
                      "d18o_03", "d18o_se_03", "d18o_04", "d18o_se_04",
                      "d18o_05", "d18o_se_05", "d18o_06", "d18o_se_06",
                      "d18o_07", "d18o_se_07", "d18o_08", "d18o_se_08",
                      "d18o_09", "d18o_se_09", "d18o_10", "d18o_se_10",
                      "d18o_11", "d18o_se_11", "d18o_12", "d18o_se_12")
           eType = 2
         },
         "USSurf" = {
           dpath.post = "USSws.zip" 
           lnames = c("d2h.tif", "d2h_se.tif", "d18o.tif", 
                      "d18o_se.tif")
           onames = c("d2h", "d2h_sd", "d18o", "d18o_sd")
           eType = 2
         },
         "USTap" = {
           dpath.post = "USTap.zip" 
           lnames = c("d2h.tif", "d2h_se.tif", "d2h_sd.tif", "d18o.tif", 
                      "d18o_se.tif", "d18o_sd.tif")
           onames = c("d2h", "d2h_se", "d2h_sd", "d18o", "d18o_se",
                      "d18o_sd")
           eType = 2
         },
         "USSr" = {
           dpath.post = "USSr.zip" 
           lnames = c("USSr_Rock.tif", "USSr_Weath.tif", 
                      "USSr_Riv.tif")
           onames = c("sr_rock", "sr_weath", "sr_riv")
           eType = 1
         },
         "CaribSr" = {
           dpath.post = "CaribSr.zip" 
           lnames = c("CaribSr_Rock.tif", "CaribSr_Weath.tif",
                      "CaribSr_Riv.tif")
           onames = c("sr_rock", "sr_weath", "sr_riv")
           eType = 1
         },
         stop("isoType invalid"))

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
  
  if(!file.exists(dpath.post)){
    dfs = dlf(paste0(dpath.pre, dpath.post), dpath.post, ot, wd)
    pdlf(dfs, wd, ot, isoType)
  }
  
  procRest = function(fn, lnames){
    if(!all(lnames %in% list.files())){
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
    procRest(dpath.post, lnames)    
  },
  error = function(cond){
    stop(cond)
  },
  finally = {
    options(timeout = ot)
    setwd(wd)
  })
  
  switch(eType,
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

