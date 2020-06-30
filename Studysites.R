library(maps)
library(mapdata)
library(raster)
library(GISTools)
library(prettymapr)

#Map of Newfoundland:

#xaxis labels
xat<-as.numeric(c("-59","-57","-55","-53"))
xlab<-(c("59˚W","57˚W","55˚ W","53 ˚W"))
ylab<-(c("46˚N","47˚N","48˚N","49˚N","50˚N","51˚N","52˚N"))

#yaxis labels one option:
yat<-as.numeric(c("48","50","52","54"))
yat <- pretty(m$range[3:4])
xat <- pretty(m$range[1:2])

## OR if you want to add the degree sign afterwards. 
#Above is if you need to have more specific with a decimal place and the W or N symbols
xlab <- parse(text=degreeLabelsEW(xat))
ylab <- parse(text=degreeLabelsNS(yat))


m<-map("worldHires","Canada", 
       xlim=c(-59.5, -51.00), 
       ylim=c(46,55), 
       col="grey75", fill=TRUE)

# to add points for a sample site
points(-54.183789, 46.82149, pch=16, col="black", cex=1) #CSM
points(-54.116889, 49.84098, pch=16, col="black", cex=1) #LFOGO
points(-56.575889, 53.93645, pch=16, col="black", cex=1) #GANN
points(-52.773818, 47.26220, pch=16, col="black", cex=1) #GULL
points(-52.809318, 47.18205, pch=16, col="black", cex=1) #GREAT

text(-54.883789, 46.62149, "Cape St. Mary's", col="black", cex=0.75) #CSM
text(-54.116889, 49.99, "Little Fogo Islands", col="black", cex=0.75) #LFOGO
text(-56.3, 54.1, "Gannet Islands", col="black", cex=0.75) #GANN
text(-51.9, 47.35220, "Gull Island", col="black", cex=0.75) #GULL
text(-51.9, 47.08005, "Great Island", col="black", cex=0.75) #GREAT

#add xaxis
axis(1, at = xat, labels=xlab,cex.axis=0.75)

#add yaxis
axis(2,  at=yat, labels=ylab, cex.axis=0.75)


#add a box around the map
box(lwd=0.75)

#add rectangle around study site
rect(-52.850,47.177,-52.750,47.265, col = c(NA,0),border = "black", lwd =2.25)

#north arrow, only works with some packages
addnortharrow(pos = "topright", padin = c(0.15, 0.15), scale = 0.25,
              lwd = 1, border = "black", cols = c("white", "black"),
              text.col = "black")

##small map
# Finding ISO3 code for Canada
library(raster)
getData('ISO3')  # Canada's code is "CAN"
# Reading in data at different levels
can0<-getData('GADM', country="CAN", level=0)

map(can0, xlim=c(-52.85, -52.75), ylim=c(47.177,47.269),col="gray75", fill=TRUE)
#plot(can0,xlim=c(-52.85, -52.75), ylim=c(47.177,47.269),col="gray75")

# to add points for a sample site
points(-52.773818, 47.26220, pch=16, col="black", cex=1) #GULL
points(-52.809318, 47.18205, pch=16, col="black", cex=1) #GREAT

text(-52.779818, 47.25220, "Gull Island", col="black", cex=1) #GULL
text(-52.80318, 47.18005, "Great Island", col="black", cex=1) #GREAT

#axis
axis(1,at = c(-52.85,-52.8,-52.75), labels=c("52.85˚W","52.8˚W","52.75˚W"))
axis(2,at = c(47.2,47.24), labels = c("47.20˚N","47.24˚N"))
#add a box around the map
box(lwd=0.75)
