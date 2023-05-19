options(stringsAsFactors = FALSE)
library(openxlsx)
library(devtools)
library(terra)
library(assignR)

#----

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
  ),
  "GlobalSr" = list(
    dpath.post = "GlobalSr.zip",
    lnames = c("GlobalSr.tif", "GlobalSr_se.tif"),
    onames = c("sr_bio", "sr_bio_se"),
    eType = 2
  )
)

#Save internal
use_data(GIconfig, internal = TRUE, overwrite = TRUE)

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

#Convert to SpatVector
knownOrig_sites = vect(sites, geom = c("Longitude", "Latitude"), crs = "WGS84")

#Write knownOrig parts
writeVector(knownOrig_sites, "inst/extdata/knownOrig_sites.shp")
write.csv(knownOrig_samples, "inst/extdata/knownOrig_samples.csv", row.names = FALSE)
write.csv(knownOrig_sources, "inst/extdata/knownOrig_sources.csv", row.names = FALSE)

stds = list(hstds = hstds, ostds = ostds, ham = ham, oam = oam)

#Write it all to /data/
use_data(stds, overwrite = TRUE)

#Prepare MI strontium isoscape
sr = getIsoscapes("USSr")
sr = sr$sr_weath
srun = setValues(sr, values(sr) * 0.01)
sr = c(sr, srun)
states.proj = project(states, crs(sr))
mi = states.proj[states.proj$STATE_NAME == "Michigan",]
sr_MI = mask(sr, mi)
sr_MI = crop(sr_MI, mi)
names(sr_MI) = c("weathered.mean", "weathered.sd")
sr_MI = aggregate(sr_MI, 10)
writeRaster(sr_MI, "inst/extdata/sr_MI.tif", overwrite = TRUE)

#Prepare lrNA H isoscape
pcp = getIsoscapes()
pcp = c(pcp$d2h, pcp$d2h.se)
pcp = mask(pcp, naMap)
pcp = crop(pcp, naMap)
d2h_lrNA = aggregate(pcp, 48, na.rm = TRUE)
crs(d2h_lrNA) = crs("WGS84")
writeRaster(d2h_lrNA, "inst/extdata/d2h_lrNA.tif", overwrite = TRUE)
