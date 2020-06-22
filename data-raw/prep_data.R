options(stringsAsFactors = FALSE)
library(openxlsx)
library(sp)
library(devtools)

#World outline map
load("data-raw/wrld_simpl.rda")

#grap WGS84 projection
p = proj4string(wrld_simpl)

#This one is internal
use_data(wrld_simpl, internal = TRUE)

#Standards adjacency matrix for H
ham = read.xlsx("data-raw/ham.xlsx", rowNames = TRUE)
ham = as.matrix(ham)
#Verify matrix symmetry
isSymmetric(ham)

#Standards adjacency matrix for O
oam = read.xlsx("data-raw/oam.xlsx", rowNames = TRUE)
oam = as.matrix(oam)
#Verify matrix xymmetry
isSymmetric(oam)

#Standards definitions files
hsds = read.xlsx("data-raw/hsds.xlsx")
osds = read.xlsx("data-raw/osds.xlsx")

#Verify rownumber matches adjacency matrix dimensions
nrow(hsds) == nrow(ham)
nrow(osds) == nrow(oam)

#Verify that all matrix entries have a match in definition file
all(row.names(ham) %in% hsds$Scale)
all(row.names(oam) %in% osds$Scale)

#Known origin data table
knownOrig_sources = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_sources")
sites = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_sites")
knownOrig_samples = read.xlsx("data-raw/knownOrigNew.xlsx", 
                              sheet = "knownOrig_samples")

#check standard scale names
ss = unique(knownOrig_sources$H_std_scale)
ss = ss[!is.na(ss)]
all(ss %in% hsds$Scale)
ss = unique(knownOrig_sources$O_std_scale)
ss = ss[!is.na(ss)]
all(ss %in% osds$Scale)

#check linking fields
all(knownOrig_samples$Site_ID %in% sites$Site_ID)
all(knownOrig_samples$Dataset_ID %in% knownOrig_sources$Dataset_ID)

#Convert to SPDF
knownOrig_sites = SpatialPointsDataFrame(sites[,2:3], 
                                         data = sites[,c(1,4:ncol(sites))],
                                         proj4string = CRS(p))
  
#Write it all to /data/
use_data(ham, oam, hsds, osds, knownOrig_samples, knownOrig_sites, 
          knownOrig_sources, overwrite = TRUE)
