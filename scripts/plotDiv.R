setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")
load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/currentDiv.Rdata")

library(maptools)
data(wrld_simpl)
# Transform World as raster
r <- raster(wrld_simpl, res = 0.5)
wrld_r <- rasterize(wrld_simpl, r)

pdf("plots/maps.pdf")

###
# First lot
###

for(i in 1:(length(list_r.div.mean)/3)){
  
#get some of all r.div.mean among indivduals
  #could also do /max?
sumRas<-sum(list_r.div.mean[[i]]@data@values,na.rm = T)

#get proportion of toatl r.div.mean in each pairwise comparison
propRas<-list_r.div.mean[[i]]@data@values/sumRas
list_r.div.mean[[i]]@data@values<-propRas

if(i == 1){
  mos<-list_r.div.mean[[i]]
} else {
  mos<-mosaic(mos,list_r.div.mean[[i]],fun='sum')
}
}

par(mar=c(4,4,2,2))
par(oma=c(0,0,0,2))
plot(wrld_r, col = "grey",xlim=c(-20,30),ylim=c(-10,15),legend=F,xlab="Lon",ylab="Lat")
mos@data@values[mos@data@values<0.00001] = NA
plot(mos,add=T)
mtext(line = 2, side = 4,"Genetic diversity")


###
# Second lot
###

for(i in ((length(list_r.div.mean)/3)+1):(length(list_r.div.mean)/3)*2){
  
  #get some of all r.div.mean among indivduals
  #could also( do /max?
  sumRas<-sum(list_r.div.mean[[i]]@data@values,na.rm = T)
  
  #get proportion of toatl r.div.mean in each pairwise comparison
  propRas<-list_r.div.mean[[i]]@data@values/sumRas
  list_r.div.mean[[i]]@data@values<-propRas
  
  if(i == 1){
    mos<-list_r.div.mean[[i]]
  } else {
    mos<-mosaic(mos,list_r.div.mean[[i]],fun='sum')
  }
}

par(mar=c(4,4,2,2))
par(oma=c(0,0,0,2))
plot(wrld_r, col = "grey",xlim=c(-20,30),ylim=c(-10,15),legend=F,xlab="Lon",ylab="Lat")
mos@data@values[mos@data@values<0.00001] = NA
plot(mos,add=T)
mtext(line = 2, side = 4,"Genetic diversity")


###
# Third lot
###

for(i in (((length(list_r.div.mean)/3)*2)+1):((length(list_r.div.mean)/3)*3)){
  
  #get some of all r.div.mean among indivduals
  #could also( do /max?
  sumRas<-sum(list_r.div.mean[[i]]@data@values,na.rm = T)
  
  #get proportion of toatl r.div.mean in each pairwise comparison
  propRas<-list_r.div.mean[[i]]@data@values/sumRas
  list_r.div.mean[[i]]@data@values<-propRas
  
  if(i == 1){
    mos<-list_r.div.mean[[i]]
  } else {
    mos<-mosaic(mos,list_r.div.mean[[i]],fun='sum')
  }
}

par(mar=c(4,4,2,2))
par(oma=c(0,0,0,2))
plot(wrld_r, col = "grey",xlim=c(-20,30),ylim=c(-10,15),legend=F,xlab="Lon",ylab="Lat")
mos@data@values[mos@data@values<0.00001] = NA
plot(mos,add=T)
mtext(line = 2, side = 4,"Genetic diversity")


dev.off()














load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/currentDiv.Rdata")

##
# Max val
##
for(i in 1:length(list_r.div.mean)){
  
  #get some of all r.div.mean among indivduals
  #could also do /max?
  sumRas<-max(list_r.div.mean[[i]]@data@values,na.rm = T)
  
  #get proportion of toatl r.div.mean in each pairwise comparison
  propRas<-list_r.div.mean[[i]]@data@values/sumRas
  list_r.div.mean[[i]]@data@values<-propRas
  
  if(i == 1){
    mos<-list_r.div.mean[[i]]
  } else {
    mos<-mosaic(mos,list_r.div.mean[[i]],fun='sum')
  }
}

plot(mos,cex=0.5)
plot('worldHiRes',add=T)


plot(wrld_r, col = "grey",xlim=c(-20,40),ylim=c(-10,15))
mos@data@values[mos@data@values<0.00001] = NA
plot(mos,add=T)
