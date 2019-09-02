rm(list=ls())
library(adegenet)
library(pegas)
#library(PopGenome)
#library(poppr)

setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")
source("localDiv.r")

###
#Geographic data
###

id_ind<-read.csv("~/Downloads/ID_index.csv")
id_vou<-read.csv("~/Downloads/ID_voucher.csv")

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
foo<-foo[order(foo$index),]

#extract geographic data
#name | lon | lat

datageo<-data.frame(foo$ddlong,foo$ddlat)
rownames(datageo)<-foo$index
datageo

colnames(datageo)<-c("lon","lat")

#calculate genetic distances
#great-circle distances (km)
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) 

###
#MICROSATS
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,-1:-2],sep="-",NA.char="0",type = "codom")

###
# SNP
###

library(pegas)
library(vcfR)
vcf <- read.vcfR("~/Downloads/anni_filtered.sub2.vcf",checkFile = T, convertNA = T) #read in all data
head(vcf) 
vcf

### convert to genlight
dataind <- vcfR2genind(vcf)

#compute genetic distances
mHij<-compGen(dataind,refDglobal = 1000)

#example of calling the function "localDiv" to compute local diversity 
localDiv(latRef=2,lonRef=13,Radius=200,refD=10,datageo=datageo,geodist=geodist,mHij=mHij)

#scan diversity
sD<-scanDiv(Dscan = 0.5, Radius = 200, refD = 10)

#plot raster + data points
library(rgdal)
library(raster)

r.div.refD=r.div.mean <- raster(xmn=min(sD$Div.table$lon)-sD$Dscan/2, xmx=max(sD$Div.table$lon)+sD$Dscan/2, 
            ymn=min(sD$Div.table$lat)-sD$Dscan/2, ymx=max(sD$Div.table$lat)+sD$Dscan/2,resolution=sD$Dscan)
ncell(r.div.refD)

#plot div for ref Dist
r.div.refD<-setValues(r.div.refD,sD$Div.table.cleaned$div.refD)
plot(r.div.refD,addfun=points(datageo$lon,datageo$lat,cex=0.2))
title(paste("div.refD"),sub=paste0("Radius=",sD$Radius,"km, refD=",sD$refD,"km, ScanResol=",sD$Dscan,"?, Ncells=",sum(!is.na(sD$Div.table.cleaned$div.refD))))

#plot mean divergence
r.div.mean<-setValues(r.div.mean,sD$Div.table$div.mean)
plot(r.div.mean,addfun=points(datageo$lon,datageo$lat,cex=0.2))
title(paste("div.mean"),sub=paste0("Radius=",sD$Radius,"km, refD=",sD$refD,"km, ScanResol=",sD$Dscan,"?, Ncells=",sum(!is.na(sD$Div.table$div.mean))))

library(maps)
library(plotrix)
library(scales)
library(raster)
library(mapdata)

cmr <- getData("GADM", country='CMR',level = 0) # level = 1 provides province data

map('worldHires',xlim = c(-20,30),ylim=c(-10,30))
plot(r.div.mean,addfun=points(datageo$lon,datageo$lat,cex=0.2))
          