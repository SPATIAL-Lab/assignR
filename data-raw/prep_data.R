options(stringsAsFactors = FALSE)
library(openxlsx)
library(sp)
library(devtools)
library(raster)
library(assignR)

#WGS84 projection
p = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

#----

#World outline map
load("data-raw/wrld_simpl.rda")

#Download configuration data
map_ln = c("d2h_MA.tif", "d2h_se_MA.tif", 
           "d18o_MA.tif", "d18o_se_MA.tif")

map_on = c("d2h_MA", "d2h_se_MA", "d18o_MA", "d18o_se_MA")

mop_ln = c("d2h_MA.tif", "d2h_se_MA.tif",
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

mop_on = c("d2h_01", "d2h_se_01", "d2h_02", "d2h_se_02", 
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

GIconfig = list(
  "GlobalPrecipGS" = list(
    dpath.post = "GlobalPrecipGS.zip",
    lnames = c("d2h_GS.tif", "d2h_se_GS.tif", 
               "d18o_GS.tif", "d18o_se_GS.tif"),
    onames = c("d2h", "d2h.se", "d18o", "d18o.se"),
    eType = 2
  ),
  "USPrecipMA" = list(
    dpath.post = "USPrecip.zip",
    lnames = map_ln,
    onames = map_on,
    eType = 2
  ),
  "GlobalPrecipMA" = list(
    dpath.post = "GlobalPrecip.zip",
    lnames = map_ln,
    onames = map_on,
    eType = 2
  ),
  "USPrecipMO" = list(
    dpath.post = "USPrecip.zip",
    lnames = mop_ln,
    onames = mop_on,
    eType = 2
  ),
  "GlobalPrecipMO" = list(
    dpath.post = "GlobalPrecip.zip",
    lnames = mop_ln,
    onames = mop_on,
    eType = 2
  ),
  "USPrecipALL" = list(
    dpath.post = "USPrecip.zip",
    lnames = c(map_ln, mop_ln),
    onames = c(map_on, mop_on),
    eType = 2
  ),
  "GlobalPrecipALL" = list(
    dpath.post = "GlobalPrecip.zip",
    lnames = c(map_ln, mop_ln),
    onames = c(map_on, mop_on),
    eType = 2
  ),
  "USSurf" = list(
    dpath.post = "USSw.zip",
    lnames = c("d2h.tif", "d2h_se.tif", "d18o.tif", 
               "d18o_se.tif"),
    onames = c("d2h", "d2h_sd", "d18o", "d18o_sd"),
    eType = 2
  ),
  "USTap" = list(
    dpath.post = "USTap.zip",
    lnames = c("d2h.tif", "d2h_se.tif", "d2h_sd.tif", "d18o.tif", 
               "d18o_se.tif", "d18o_sd.tif"),
    onames = c("d2h", "d2h_se", "d2h_sd", "d18o", "d18o_se",
               "d18o_sd"),
    eType = 2
  ),
  "USSr" = list(
    dpath.post = "USSr.zip",
    lnames = c("USSr_Rock.tif", "USSr_Weath.tif", 
               "USSr_Riv.tif"),
    onames = c("sr_rock", "sr_weath", "sr_riv"),
    eType = 1
  ),
  "CaribSr" = list(
    dpath.post = "CaribSr.zip",
    lnames = c("CaribSr_Rock.tif", "CaribSr_Weath.tif",
               "CaribSr_Riv.tif"),
    onames = c("sr_rock", "sr_weath", "sr_riv"),
    eType = 1
  )
)

#Save internal
use_data(wrld_simpl, GIconfig, internal = TRUE, overwrite = TRUE)

#adjacency matrix for H
ham = read.xlsx("data-raw/ham.xlsx", rowNames = TRUE)
ham = as.matrix(ham)
#Verify matrix symmetry
isSymmetric(ham)

#adjacency matrix for O
oam = read.xlsx("data-raw/oam.xlsx", rowNames = TRUE)
oam = as.matrix(oam)
#Verify matrix xymmetry
isSymmetric(oam)

#Standards definitions files
hstds = read.xlsx("data-raw/hstds.xlsx")
ostds = read.xlsx("data-raw/ostds.xlsx")

#Verify rownumber matches adjacency matrix dimensions
nrow(hstds) == nrow(ham)
nrow(ostds) == nrow(oam)

#Verify that all matrix entries have a match in definition file
all(row.names(ham) %in% hstds$Calibration)
all(row.names(oam) %in% ostds$Calibration)

#Known origin data table
knownOrig_sources = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_sources")
sites = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_sites")
knownOrig_samples = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_samples")

#check standard scale names
ss = unique(knownOrig_sources$H_cal)
ss = ss[!is.na(ss)]
all(ss %in% hstds$Calibration)
ss = unique(knownOrig_sources$O_cal)
ss = ss[!is.na(ss)]
all(ss %in% ostds$Calibration)

#check linking fields
all(knownOrig_samples$Site_ID %in% sites$Site_ID)
all(knownOrig_samples$Dataset_ID %in% knownOrig_sources$Dataset_ID)

#Convert to SPDF
knownOrig_sites = SpatialPointsDataFrame(sites[,2:3], 
                                         data = sites[,c(1,4:ncol(sites))],
                                         proj4string = p)

#update to include WKT representation; requires rgdal >=1.5-17
knownOrig_sites = rebuild_CRS(knownOrig_sites)

#Group data objects
knownOrig = list(sites = knownOrig_sites, samples = knownOrig_samples, 
                 sources = knownOrig_sources)

stds = list(hstds = hstds, ostds = ostds, ham = ham, oam = oam)
  
#Prepare MI strontium isoscape
sr = getIsoscapes("USSr")
sr = sr$sr_weath
srun = setValues(sr, getValues(sr) * 0.01)
sr = brick(sr, srun)
states.proj = spTransform(states, crs(sr))
mi = states.proj[states.proj$STATE_NAME == "Michigan",]
sr_MI = mask(sr, mi)
sr_MI = crop(sr_MI, mi)
names(sr_MI) = c("weathered.mean", "weathered.sd")
sr_MI = aggregate(sr_MI, 10)

#Write it all to /data/
use_data(knownOrig, stds, sr_MI, overwrite = TRUE)
