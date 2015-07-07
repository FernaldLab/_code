#using pitsmall
library(MBA);
library(rgl);
#load('~/Downloads/Bower Plots.RData')
pitcastle=as.matrix(read.table('~/Downloads/pitcastlecleaned.asc'))

#outk=kmeans(pitsmall,nrow(pitsmall)/10);
outk=kmeans(pitcastle,nrow(pitcastle)/1000);

reduced=outk$centers

outreduced=mba.surf(reduced,nrow(reduced),nrow(reduced))

tmp=outreduced$xyz.est


persp3d(tmp, theta = 135, phi = 120, col = "green3", scale = FALSE,
  ltheta = -120, shade = 1.5, expand = 1, axes=T, box = T)