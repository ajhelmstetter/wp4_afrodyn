setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")
load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/futureDiv.Rdata")

list_r.div.mean_fut<-list_r.div.mean

load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/currentDiv.Rdata")

names(list_r.div.mean)
names(list_r.div.mean_fut)

plot(list_r.div.mean[[10]])
plot(list_r.div.mean_fut[[10]])

list_r.div.mean<-list_r.div.mean[names(list_r.div.mean)%in%names(list_r.div.mean_fut)]
names(list_r.div.mean)
names(list_r.div.mean_fut)


diff<-re-list_r.div.mean[[2]]
plot(diff)

e<-extent(list_r.div.mean[[2]])
re <- extend(list_r.div.mean_fut[[2]], e)


pdf("plots/gen_diff.pdf")

for(i in 1:(length(list_r.div.mean)/3)){

list_r.div.mean[[i]][is.na(list_r.div.mean[[i]])] <- 0 
list_r.div.mean_fut[[i]][is.na(list_r.div.mean_fut[[i]])] <- 0 

e<-extent(list_r.div.mean[[i]])
re <- extend(list_r.div.mean_fut[[i]], e)

diff<-re-list_r.div.mean[[i]]

library(rasterVis)
par(mar=c(1,1,1,1))
print(levelplot(diff,par.settings = RdBuTheme))

print(levelplot(diff, layers = 1, contour=TRUE))

#vectorplot(diff, par.settings=RdBuTheme())

}

dev.off()

library(maptools)
data(wrld_simpl)
# Transform World as raster
r <- raster(wrld_simpl, res = 0.5)
wrld_r <- rasterize(wrld_simpl, r)

pdf("plots/future_maps.pdf")

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

###
# 
###

str(list_gen_dist)

curr_mean_dist<-vector()
curr_sum_dist<-vector()

for(i in 1:length(list_gen_dist)){
  curr_mean_dist[i]<-mean(list_gen_dist[[i]])
  curr_sum_dist[i]<-sum(list_gen_dist[[i]])
}

names(curr_mean_dist)<-names(list_gen_dist)
names(curr_sum_dist)<-names(list_gen_dist)

curr_mean_dist
curr_sum_dist


fut_mean_dist<-vector()
fut_sum_dist<-vector()

for(i in 1:length(list_gen_dist_fut)){
  fut_mean_dist[i]<-mean(list_gen_dist_fut[[i]])
  fut_sum_dist[i]<-sum(list_gen_dist_fut[[i]])
}

names(fut_mean_dist)<-names(list_gen_dist_fut)
names(fut_sum_dist)<-names(list_gen_dist_fut)

fut_mean_dist
fut_sum_dist











load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/futureDiv.Rdata")

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
