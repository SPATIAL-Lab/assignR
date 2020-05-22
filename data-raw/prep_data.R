options(stringsAsFactors = FALSE)

#World outline map
load("data-raw/wrld_simpl.rda")

#This one is internal
use_data(wrld_simpl, internal = TRUE)

#Standards adjacency matrix for H
ham = read.csv("data-raw/H_connections.csv", row.names = 1)
ham = as.matrix(ham)
#Verify matrix symmetry
isSymmetric(ham)

#Standards adjacency matrix for O
oam = read.csv("data-raw/O_connections.csv", row.names = 1)
oam = as.matrix(oam)
#Verify matrix xymmetry
isSymmetric(oam)

#Standards definitions files
hsds = read.csv("data-raw/H_stds.csv")
hsds = hsds[, 1:17]
osds = read.csv("data-raw/O_stds.csv")
osds = osds[, 1:17]
#Verify rownumber matches adjacency matrix dimensions
nrow(hsds) == nrow(ham)
nrow(osds) == nrow(oam)
#Verify that all matrix entries have a match in definition file
all(row.names(ham) %in% hsds$Scale)
all(row.names(oam) %in% osds$Scale)

#Known origin data table
knownOrig
#Convert to SPDF

  
#Known origin projects table
projects 
  
#Write it all to /data/
use_data(ham, oam, hsds, osds, overwrite = TRUE)
