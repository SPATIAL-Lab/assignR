
### Data transformation through assignR function subOrigData & analysis of pre- and post-transformation data 

## Install assignR development release
library(devtools)
#install_github("SPATIAL-Lab/assignR@*release") #will be bundled as a release when we are finished and ready to resubmit
install_github("SPATIAL-Lab/assignR")

## Load assignR & other required libraries
library(assignR)
library(raster)
library(sp)

## Download precip isoscapes
download.file("http://wateriso.utah.edu/waterisotopes/media/ArcGrids/IsotopeMaps.zip",
                        "pcp_iso.zip")
unzip("pcp_iso.zip")

## Load mean annual d2H and d18O
## Check against file paths output above which may differ on your system 
H = raster("./Hma.asc", gz=FALSE)
O = raster("./Oma.asc", gz=FALSE)

## Load the knownOrig dataset, which is a list of three data objects: sources (data frame), sites (spatial point data frame),
# and samples (data frame)
data("knownOrig")

#class(knownOrig)
#str(knownOrig)

#class(knownOrig$sources) 
#class(knownOrig$sites) 
#class(knownOrig$samples)

#names(knownOrig$sources)
#names(knownOrig$sites)
#names(knownOrig$samples)


## Boxplots of original and recalibrated data 

## From the knownOrig dataset subset datasets analyzing keratins using comparative equilibration methods  
DatasetIDs_d2H = c(9,10,8,19,11,13,5,14,2,6,17,12,7,15,16)
DatasetIDs_d18O = c(9,10,8,5,7,15,16)

## Extract untransformed data calibrated to original reference scales
Origd2H = subOrigData(marker = "d2H", dataset = DatasetIDs_d2H, ref_scale = NULL)  
Recald2H = subOrigData(marker = "d2H", dataset = DatasetIDs_d2H, ref_scale = "VSMOW_H") 

## Recalibrate data to the VSMOW reference scale and extract transformed data  
Origd18O = subOrigData(marker = "d18O", dataset = DatasetIDs_d18O, ref_scale = NULL)
Recald18O = subOrigData(marker = "d18O", dataset = DatasetIDs_d18O, ref_scale = "VSMOW_O") 


## Bind original and recalibrated data by row and define their category (Orig vs. Recal) 
Origd2H = Origd2H$data
Recald2H = Recald2H$data
Origd18O = Origd18O$data
Recald18O = Recald18O$data

Origd2H$Cat = rep("Orig",length(Origd2H$d2H))
Recald2H$Cat = rep("Recal",length(Recald2H$d2H))
Origd18O$Cat = rep("Orig",length(Origd18O$d18O))
Recald18O$Cat = rep("Recal",length(Recald18O$d18O))

d2H = rbind(Origd2H, Recald2H)
d18O = rbind(Origd18O, Recald18O)


## Set order in which the datasets should be plotted
d2H$Dataset_ID = factor(d2H$Dataset_ID, levels=c(9,10,8,19,11,13,5,14,2,6,17,12,7,15,16))
d18O$Dataset_ID = factor(d18O$Dataset_ID, levels=c(9,10,8,5,7,15,16))

## Set colors for plotting
cols = gray.colors(n=6, gamma = 1)

## Fig. 1A 
pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.1A_final.pdf", height=5, width=8, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.1A_final.pdf", height=5, width=8, encoding="WinAnsi.enc")

par(mfrow=c(1,1), oma=c(0,0,0.5,0.5), mar=c(7.5,4,0.5,0.5), mgp=c(2.5,1,0)) 

boxplot(d2H ~ Cat:Dataset_ID, data= d2H, col=c(cols[2],cols[5]), ylim=c(-200,50), xaxt="n", ylab=expression(paste("Keratin ", delta^{2}, "H (\u2030)")), xlab="", cex.axis = 1, cex.lab=1, outline=FALSE)
labels=paste(c("Bowen et al. 2009","Ehleringer et al. 2008","Thompson et al. 2010","Bataille Human Hair","Wunder Plover","Neto et al. 2006","Hobson et al. 2004","Hobson et al. 2009","Hobson et al. 2012","Lott & Smith 2006","Prochazka et al. 2013","van Dijk et al. 2014","Hobson & Kohler 2015","Bowen Scaup","Magozzi Towhee"))
axis(side=1, at=1.5+2*(0:14), labels=FALSE)
box(lwd=1.25)
text(x = 1.5+2*(0:14), y = par("usr")[3]-9, srt=60, adj=1, xpd=TRUE, labels=labels, cex=0.75)

points(15, -108, pch=16, col="black", cex=1)
points(16, -80.82465, pch=16, col="black", cex=1)
points(27, -110, pch=16, col="black", cex=1)
points(28, -79.66745, pch=16, col="black", cex=1)

legend("topright", c("Original","Recalibrated"), fill=c(cols[2], cols[5]), bty = "n", cex=0.75)

dev.off()

## Fig. 1B
pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.1B_final.pdf", height=5, width=8, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.1B_final.pdf", height=5, width=8, encoding="WinAnsi.enc")

par(mfrow=c(1,1), oma=c(0,0,0.5,0.5), mar=c(7.5,4,0.5,0.5), mgp=c(2.5,1,0)) 

boxplot(d2H.sd ~ Cat:Dataset_ID, data= d2H, col=c(cols[2],cols[5]), ylim=c(0,5), xaxt="n", ylab=expression(paste("Keratin SE ", delta^{2}, "H (\u2030)")), xlab="", cex.axis = 1, cex.lab=1, outline=FALSE)
labels=paste(c("Bowen et al. 2009","Ehleringer et al. 2008","Thompson et al. 2010","Bataille Human Hair","Wunder Plover","Neto et al. 2006","Hobson et al. 2004","Hobson et al. 2009","Hobson et al. 2012","Lott & Smith 2006","Prochazka et al. 2013","van Dijk et al. 2014","Hobson & Kohler 2015","Bowen Scaup","Magozzi Towhee"))
axis(side=1, at=1.5+2*(0:14), labels=FALSE)
box(lwd=1.25)
text(x = 1.5+2*(0:14), y = par("usr")[3]-0.20, srt=60, adj=1, xpd=TRUE, labels=labels, cex=0.75)

points(15, 0, pch=16, col="black", cex=1)
points(16, 0.3553910, pch=16, col="black", cex=1)
points(27, 0, pch=16, col="black", cex=1)
points(28, 0.6865116, pch=16, col="black", cex=1)

dev.off()

## Fig. 2A
pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.2A_final.pdf", height=5, width=8, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.2A_final.pdf", height=5, width=8, encoding="WinAnsi.enc")

par(mfrow=c(1,1), oma=c(0,0,0.5,0.5), mar=c(7.5,4,0.5,0.5), mgp=c(2.5,1,0)) 

boxplot(d18O ~ Cat:Dataset_ID, data= d18O, col=c(cols[2],cols[5]), ylim=c(-5,23), 
        xaxt="n", ylab=expression(paste("Keratin ", delta^{18}, "O (\u2030)")), xlab="", cex.axis = 1, cex.lab=1, outline=FALSE)
labels=paste(c("Bowen et al. 2009","Ehleringer et al. 2008","Thompson et al. 2010","Hobson et al. 2004","Hobson & Kohler 2015","Bowen Scaup","Magozzi Towhee"))
axis(side=1, at=1.5+2*(0:6), labels=FALSE)
box(lwd=1.25)
text(x = 1.5+2*(0:6), y = par("usr")[3]-1, srt=60, adj=1, xpd=TRUE, labels=labels, cex=0.75)

dev.off()

## Fig 2B
pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.2B_final.pdf", height=5, width=8, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.2B_final.pdf", height=5, width=8, encoding="WinAnsi.enc")

par(mfrow=c(1,1), oma=c(0,0,0.5,0.5), mar=c(7.5,4,0.5,0.5), mgp=c(2.5,1,0)) 

boxplot(d18O.sd ~ Cat:Dataset_ID, data= d18O, col=c(cols[2],cols[5]), ylim=c(0,0.5),  
        xaxt="n", ylab=expression(paste("Keratin SE ", delta^{18}, "O (\u2030)")), xlab="", cex.axis = 1, cex.lab=1, outline=FALSE)
labels=paste(c("Bowen et al. 2009","Ehleringer et al. 2008","Thompson et al. 2010","Hobson et al. 2004","Hobson & Kohler 2015","Bowen Scaup","Magozzi Towhee"))
axis(side=1, at=1.5+2*(0:6), labels=FALSE)
box(lwd=1.25)
text(x = 1.5+2*(0:6), y = par("usr")[3]-0.025, srt=60, adj=1, xpd=TRUE, labels=labels, cex=0.75)

dev.off()


## Tissue-water relationships with original and recalibrated data  

## Bind original and recalibrated data by columns
d2H.lm = Origd2H
d2H.lm$Recald2H = Recald2H$d2H
d2H.lm$Recald2H.sd = Recald2H$d2H.sd

## Match Dataset IDs with Dataset names
d2H.lm$Dataset_name = ifelse(d2H.lm$Dataset_ID == 2, "Hobson et al. 2012 Plos",
                             ifelse(d2H.lm$Dataset_ID == 7, "Hobson and Kohler 2015 Ecol Evol",
                                    ifelse(d2H.lm$Dataset_ID == 5, "Hobson et al. 2004 Oecologia",      
                                           ifelse(d2H.lm$Dataset_ID == 6, "Lott and Smith 2006 Auk",       
                                                  ifelse(d2H.lm$Dataset_ID == 8, "Thompson et al. 2010 Am J Phys Anthropol",
                                                         ifelse(d2H.lm$Dataset_ID == 9, "Bowen et al. 2009 Am J Phys Anthropol",       
                                                                ifelse(d2H.lm$Dataset_ID == 10, "Ehleringer et al. 2008 PNAS",
                                                                       ifelse(d2H.lm$Dataset_ID == 11, "Wunder Plover",
                                                                              ifelse(d2H.lm$Dataset_ID == 12, "van Dijk et al. 2014 J Avian Biol",
                                                                                     ifelse(d2H.lm$Dataset_ID == 13, "Neto et al. 2006 J Avian Biol",
                                                                                            ifelse(d2H.lm$Dataset_ID == 14, "Hobson et al. 2009 Plos", 
                                                                                                   ifelse(d2H.lm$Dataset_ID == 15, "Bowen Scaup",        
                                                                                                          ifelse(d2H.lm$Dataset_ID == 16, "Magozzi Towhee",       
                                                                                                                 ifelse(d2H.lm$Dataset_ID == 17, "Prochazka et al. 2013 J Avian Biol",
                                                                                                                        "Bataille Canadian Human Hair"
                       ))))))))))))))

## Remove datasets that only have state-geographic info 
d2H.lm = d2H.lm[!d2H.lm$Dataset_ID %in% c(14,15),]

## Calculate site-averages
Sites2 = as.data.frame(knownOrig$sites)

Lon = 0
Lat = 0
Site_ID = 0

for (i in 1:length(d2H.lm$Site_ID)) {
  
  Site_ID[i] = d2H.lm$Site_ID[i]
  
  Lon[i] = Sites2$Longitude[Sites2$Site_ID == Site_ID[i]]
  Lat[i] = Sites2$Latitude[Sites2$Site_ID == Site_ID[i]]
  
}

d2H.lm$Lon = Lon
d2H.lm$Lat = Lat


d2H.lm = d2H.lm[,-c(1,2,3,4,5,6,7,8,11,12,14,17,18,19,20)]
d2H.lm$Lon = round(d2H.lm$Lon, digits=3)
d2H.lm$Lat = round(d2H.lm$Lat, digits=3)

d2H.lm.avg = aggregate(.~Dataset_name+Taxon+Group+Material_type+Lon+Lat, data=d2H.lm, mean)


## Extract ann avg precip d2H values at each site and add it to the sample data frame
points = d2H.lm.avg[,5:6]
coordinates(points) = ~ Lon + Lat
projection(H) =  CRS("+proj=longlat +datum=WGS84")
d2Hpoints = extract(H, points)
pointdata = as.data.frame(cbind(coordinates(points), d2Hpoints))
colnames(pointdata) = c("Longitude", "Latitude", "d2H_precip")

d2H.lm.avg$d2H_precip = pointdata$d2H_precip


## Repeat for O 
d18O.lm = Origd18O
d18O.lm$Recald18O = Recald18O$d18O
d18O.lm$Recald18O.sd = Recald18O$d18O.sd

d18O.lm$Dataset_name = ifelse(d18O.lm$Dataset_ID == 7, "Hobson and Kohler 2015 Ecol Evol",
                             ifelse(d18O.lm$Dataset_ID == 5, "Hobson et al. 2004 Oecologia",
                                    ifelse(d18O.lm$Dataset_ID == 8, "Thompson et al. 2010 Am J Phys Anthropol",      
                                           ifelse(d18O.lm$Dataset_ID == 9, "Bowen et al. 2009 Am J Phys Anthropol",       
                                                  ifelse(d18O.lm$Dataset_ID == 10, "Ehleringer et al. 2008 PNAS",
                                                         ifelse(d18O.lm$Dataset_ID == 15, "Bowen Scaup",       
                                                                "Magozzi Towhee"
                          ))))))

d18O.lm = d18O.lm[!d18O.lm$Dataset_ID == 15,]


Lon = 0
Lat = 0
Site_ID = 0

for (i in 1:length(d18O.lm$Site_ID)) {
  
  Site_ID[i] = d18O.lm$Site_ID[i]
  
  Lon[i] = Sites2$Longitude[Sites2$Site_ID == Site_ID[i]]
  Lat[i] = Sites2$Latitude[Sites2$Site_ID == Site_ID[i]]
  
}

d18O.lm$Lon = Lon
d18O.lm$Lat = Lat


d18O.lm = d18O.lm[,-c(1,2,3,4,5,6,7,8,11,12,14,15,16,19,20)]
d18O.lm$Lon = round(d18O.lm$Lon, digits=3)
d18O.lm$Lat = round(d18O.lm$Lat, digits=3)

d18O.lm.avg = aggregate(.~Dataset_name+Taxon+Group+Material_type+Lon+Lat, data=d18O.lm, mean)


points = d18O.lm.avg[,5:6]
coordinates(points) = ~ Lon + Lat
projection(O) =  CRS("+proj=longlat +datum=WGS84")
d18Opoints = extract(O, points)
pointdata = as.data.frame(cbind(coordinates(points), d18Opoints))
colnames(pointdata) = c("Longitude", "Latitude", "d18O_precip")

d18O.lm.avg$d18O_precip = pointdata$d18O_precip



## Set plotting parameters and plot d2H relationships 
d2H.lm.avg$Group = factor(d2H.lm.avg$Group)
d2H.lm.avg = d2H.lm.avg[d2H.lm.avg$Group %in% c("Modern human","Passerine","Ground bird","Water bird", "Raptor"),]
d2H.lm.avg$Dataset_name = factor(d2H.lm.avg$Dataset_name)
dataset_names = levels(d2H.lm.avg$Dataset_name)
dataset_symbols = array(data=c(0:length(dataset_names)), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))
dataset_colors = array(data=c("darkmagenta","chartreuse3","firebrick",rep("tomato2",3),"springgreen4","orange red","tomato2","chartreuse3","tomato2","orange"), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))
dataset_size = array(data=c(rep(1,6),3,rep(1,5)), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))

pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.3_final.pdf", height=12, width=10, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.3_final.pdf", height=12, width=10, encoding="WinAnsi.enc")

par(mfrow=c(5,3), oma=c(3,3,0.5,0.5), mar=c(2,2,0.5,0.5), mgp=c(2.5,1,0)) 


## Modern human 

Mod_human = d2H.lm.avg[d2H.lm.avg$Group == "Modern human",]

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), ylab="", xaxt='n', cex.axis = 1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Mod_human$d2H_precip[Mod_human$Dataset_name == i], Mod_human$d2H[Mod_human$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Mod_human$d2H ~ Mod_human$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Mod_human$d2H_precip[Mod_human$Dataset_name == i], Mod_human$Recald2H[Mod_human$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Mod_human$Recald2H ~ Mod_human$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), xaxt='n', yaxt='n', ann=FALSE, bty='n')

## legend scales 
legend("topleft", legend=c("OldUT_H_1","CAN_H_1","DEN_H_1","OldSA.2_H_1","OldSA.1_H_1","SA_H_5","UT_H_2"), fill = c("chartreuse3","darkmagenta","orange","orange red","tomato2","firebrick","springgreen4"), bty="n", cex=1)

## legend datasets 
legend("topright", legend=c("Bataille Human Hair","Ehleringer et al. 2008","Hobson & Kohler 2015","Hobson et al. 2004","Hobson et al. 2012","Lott & Smith 2006","Magozzi Towhee","Neto et al. 2006","Prochazka et al. 2013","Thompson et al. 2010","van Dijk et al. 2014","Wunder Plover"), pch = dataset_symbols, bty="n", cex=1)


## Passerine 

Pass = d2H.lm.avg[d2H.lm.avg$Group == "Passerine",]

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Pass$d2H, Pass$Recald2H),na.rm=T), max(c(Pass$d2H, Pass$Recald2H),na.rm=T)+20), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$d2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

i = "Hobson and Kohler 2015 Ecol Evol"	
points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$d2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])

i = "Magozzi Towhee"	
points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$d2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=3, lwd=2, col=dataset_colors[dataset_names == i])

M.old = lm(Pass$d2H ~ Pass$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Pass$d2H, Pass$Recald2H),na.rm=T), max(c(Pass$d2H, Pass$Recald2H),na.rm=T)+20), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$Recald2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

i = "Hobson and Kohler 2015 Ecol Evol"	
points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$Recald2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])

i = "Magozzi Towhee"	
points(Pass$d2H_precip[Pass$Dataset_name == i], Pass$Recald2H[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=3, lwd=2, col=dataset_colors[dataset_names == i])

M.old = lm(Pass$Recald2H ~ Pass$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), xaxt='n', yaxt='n', ann=FALSE, bty='n')


## Ground bird

Ground = d2H.lm.avg[d2H.lm.avg$Group == "Ground bird",]

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Ground$d2H, Ground$Recald2H),na.rm=T), max(c(Ground$d2H, Ground$Recald2H),na.rm=T)+20), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Ground$d2H_precip[Ground$Dataset_name == i], Ground$d2H[Ground$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Ground$d2H ~ Ground$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Ground$d2H, Ground$Recald2H),na.rm=T), max(c(Ground$d2H, Ground$Recald2H),na.rm=T)+20), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Ground$d2H_precip[Ground$Dataset_name == i], Ground$Recald2H[Ground$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Ground$Recald2H ~ Ground$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), xaxt='n', yaxt='n', ann=FALSE, bty='n')


## Water bird 

Water = d2H.lm.avg[d2H.lm.avg$Group == "Water bird",]

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Water$d2H, Water$Recald2H),na.rm=T), max(c(Water$d2H, Water$Recald2H),na.rm=T)+10), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Water$d2H_precip[Water$Dataset_name == i], Water$d2H[Water$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Water$d2H ~ Water$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Water$d2H, Water$Recald2H),na.rm=T), max(c(Water$d2H, Water$Recald2H),na.rm=T)+10), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-200,0,50), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Water$d2H_precip[Water$Dataset_name == i], Water$Recald2H[Water$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Water$Recald2H ~ Water$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), xaxt='n', yaxt='n', ann=FALSE, bty='n')


## Raptor

Raptor = d2H.lm.avg[d2H.lm.avg$Group == "Raptor",]

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Raptor$d2H, Raptor$Recald2H),na.rm=T), max(c(Raptor$d2H, Raptor$Recald2H),na.rm=T)+25), ylab="", cex.axis=1.1)

for (i in dataset_names) {
  points(Raptor$d2H_precip[Raptor$Dataset_name == i], Raptor$d2H[Raptor$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Raptor$d2H ~ Raptor$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Raptor$d2H, Raptor$Recald2H),na.rm=T), max(c(Raptor$d2H, Raptor$Recald2H),na.rm=T)+25), ylab="", cex.axis=1.1)

for (i in dataset_names) {
  points(Raptor$d2H_precip[Raptor$Dataset_name == i], Raptor$Recald2H[Raptor$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Raptor$Recald2H ~ Raptor$d2H_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-200,0), ylim=c(min(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T), max(c(Mod_human$d2H, Mod_human$Recald2H),na.rm=T)+10), xaxt='n', yaxt='n', ann=FALSE, bty='n')


mtext(text=expression(paste("Original keratin ", delta^{2}, "H (\u2030)")), side=2, line=1, outer=TRUE)
mtext(text=expression(paste("Local precipitation ", delta^{2}, "H (\u2030)")), side=1, line=1, outer=TRUE) 

dev.off()


## Plot d18O relationships 
d18O.lm.avg$Group = factor(d18O.lm.avg$Group)
d18O.lm.avg = d18O.lm.avg[d18O.lm.avg$Group %in% c("Modern human","Passerine","Ground bird"),]
d18O.lm.avg$Dataset_name = factor(d18O.lm.avg$Dataset_name)
dataset_names = levels(d18O.lm.avg$Dataset_name) 
dataset_symbols = array(data=c(1,2,3,6,9), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))
dataset_colors = array(data=c("chartreuse3","firebrick","grey","springgreen4","chartreuse3"), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))
dataset_size = array(data=c(rep(1,3),3,1), dim=c(1,length(dataset_names)), dimnames=list(NULL,dataset_names))


pdf("/Users/Sarah/Dropbox/Magozzi et al. MEE_code/Figures/Fig.4_final.pdf", height=12, width=10, encoding="WinAnsi.enc")
#pdf("/Users/Sarah/Dropbox/HO_database/paper/Fig.4_final.pdf", height=12, width=10, encoding="WinAnsi.enc")

par(mfrow=c(5,3), oma=c(3,3,0.5,0.5), mar=c(2,2,0.5,0.5), mgp=c(2.5,1,0)) 


## Modern human 

Mod_human = d18O.lm.avg[d18O.lm.avg$Group == "Modern human",]

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T), max(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T)+1), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Mod_human$d18O_precip[Mod_human$Dataset_name == i], Mod_human$d18O[Mod_human$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Mod_human$d18O ~ Mod_human$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T), max(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T)+1), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Mod_human$d18O_precip[Mod_human$Dataset_name == i], Mod_human$Recald18O[Mod_human$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Mod_human$Recald18O ~ Mod_human$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T), max(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T)+1), xaxt='n', yaxt='n', ann=FALSE, bty='n', cex.axis=1.1)


## legend scales 
legend("topleft", legend=c("OldUT_O_1","IAEA_O_1","SA_O_10","UT_O_2"), fill = c("chartreuse3","grey","firebrick","springgreen4"), bty="n", cex=1)

## legend datasets 
legend("topright", legend=c("Ehleringer et al. 2008","Hobson & Kohler 2015","Hobson et al. 2004","Magozzi Towhee","Thompson et al. 2010"), pch = dataset_symbols, bty="n", cex=1)


## Passerine 

Pass = d18O.lm.avg[d18O.lm.avg$Group == "Passerine",]

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Pass$d18O, Pass$Recald18O),na.rm=T), max(c(Pass$d18O, Pass$Recald18O),na.rm=T)+1), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$d18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

i = "Hobson and Kohler 2015 Ecol Evol"	
points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$d18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])

i = "Magozzi Towhee"	
points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$d18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=3, lwd=2, col=dataset_colors[dataset_names == i])

M.old = lm(Pass$d18O ~ Pass$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Pass$d18O, Pass$Recald18O),na.rm=T), max(c(Pass$d18O, Pass$Recald18O),na.rm=T)+1), ylab="", xaxt='n', cex.axis=1.1)
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$Recald18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

i = "Hobson and Kohler 2015 Ecol Evol"	
points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$Recald18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])

i = "Magozzi Towhee"	
points(Pass$d18O_precip[Pass$Dataset_name == i], Pass$Recald18O[Pass$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=3, lwd=2, col=dataset_colors[dataset_names == i])

M.old = lm(Pass$Recald18O ~ Pass$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T), max(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T)+1), xaxt='n', yaxt='n', ann=FALSE, bty='n', cex.axis=1.1)


## Ground bird

Ground = d18O.lm.avg[d18O.lm.avg$Group == "Ground bird",]

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Ground$d18O, Ground$Recald18O),na.rm=T), max(c(Ground$d18O, Ground$Recald18O),na.rm=T)+1), ylab="")
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Ground$d18O_precip[Ground$Dataset_name == i], Ground$d18O[Ground$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Ground$d18O ~ Ground$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)


plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Ground$d18O, Ground$Recald18O),na.rm=T), max(c(Ground$d18O, Ground$Recald18O),na.rm=T)+1), ylab="")
axis(side = 1, at = seq(-20,-2,5), labels = FALSE, tck = -0.03, cex.axis=1.1)

for (i in dataset_names) {
  points(Ground$d18O_precip[Ground$Dataset_name == i], Ground$Recald18O[Ground$Dataset_name == i], pch=dataset_symbols[dataset_names == i], cex=dataset_size[dataset_names == i], col=dataset_colors[dataset_names == i])
}

M.old = lm(Ground$Recald18O ~ Ground$d18O_precip)
abline(M.old, col="black")
cf.old = round(coef(M.old),2)
cf.old[1] = ifelse(cf.old[1] <0 , cf.old[1], paste("+", cf.old[1]))
r2.old = round(summary(M.old)$adj.r.squared, 2)
eq.old = bquote(paste("y = ", .(cf.old[2]), "x ", .(cf.old[1]), ", R"^2," = ", .(r2.old), sep=""))
mtext(eq.old,3, line=-1.5, adj=0.15, col="black", cex=0.75)

plot(x=NULL, y = NULL, xlim=c(-20,-2), ylim=c(min(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T), max(c(Mod_human$d18O, Mod_human$Recald18O),na.rm=T)+1), xaxt='n', yaxt='n', ann=FALSE, bty='n')


mtext(text=expression(paste("Original keratin ", delta^{18}, "O (\u2030)")), side=2, line=1, outer=TRUE)
mtext(text=expression(paste("Local precipitation ", delta^{18}, "O (\u2030)")), side=1, line=1, outer=TRUE)

dev.off()


