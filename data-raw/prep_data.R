options(stringsAsFactors = FALSE)
library(openxlsx)
library(sp)
library(devtools)
library(raster)
library(assignR)

#WGS84 projection
p = CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

#World outline map
load("data-raw/wrld_simpl.rda")

#This one is internal
use_data(wrld_simpl, internal = TRUE, overwrite = TRUE)

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
