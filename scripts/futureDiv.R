rm(list=ls())
library(adegenet)
library(pegas)
library(vcfR)
library(rgdal)
library(raster)

setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")
source("scripts/localDiv.R")

#load image with all data
load("~/Dropbox/projects/AJH_WP4_AFRODYN/data/datageo_fut.Rdata")

list_geo<-list_geo_fut
names(list_datageo_fut)<-names(list_datageo)
list_datageo<-list_datageo_fut
list_gen<-list_gen_fut
list_gen_dist<-list()

par(mfrow=c(4,4))
#check location of leftover points
for (i in c(1,2,4:length(list_datageo))){
  plot(list_datageo[[i]]$lon~list_datageo[[i]]$lat)
}


list_r.div.refD<-list()
list_r.div.mean<-list()
list_sD<-list()

#loop through parameters
#change from 1:number of parameters you want to loop through 
for(j in 1:3){
  
  ###
  # Define parameters 
  ### 
  
  #estimate diversity in a circular window of given Radius (in km)
  
  vRadius<-c(100,200,500)
  Radius = vRadius[j]
  
  #accounting for isolation by distance by estimating diversity as the expected divergence 
  #between individuals at a given reference distance (default refD = Radius/2)
  vrefD<-c(50,100,250)
  refD = vrefD[j]
  
  #Not sure what this is, ask Olivier - maybe just a global object for running - should be same as refD?
  vrefDglobal<-c(50,100,250)
  refDglobal = vrefDglobal[j]
  
  #centered on latRef, lonRef (in degrees decimal)
  latRef = 2
  lonRef = 13
  
  #Scan local diversity around each node of a lon/lat grid of resolution Dscan in degrees
  #0.25,0.5,1.0
  Dscan = 0.5
  
  ###
  # Prepare lists to receive results
  ### 
  
  pdfname<-paste("future","Dscan",Dscan,
                 "Radius",Radius,
                 "refDGlobal",refDglobal,
                 "refD",refD,
                 "latRef",latRef,
                 "lonRef",lonRef,".pdf",sep = "_")
  pdfname
  
  pdf(paste("plots/",pdfname,sep=""))
  
  for(i in c(1,4,5,6,7,8,9,10,11,12,13,17,18,19)){
    
    #compute genetic distances
    mHij<-compGen(list_gen[[i]],refDglobal = refDglobal)
    
    #example of calling the function "localDiv" to compute local diversity 
    localDiv(latRef=latRef,lonRef=lonRef,Radius=Radius,refD=refD,datageo=list_datageo[[i]],geodist=list_geo[[i]],mHij=mHij)
    
    #scan diversity
    #Low values for Dscan are not working well
    sD<-try(scanDiv(Dscan = 0.5, Radius = Radius, refD = refD,datageo=list_datageo[[i]],geodist=list_geo[[i]]))
    
    #plot raster + data points
    r.div.refD=r.div.mean <- raster(xmn=min(sD$Div.table$lon)-sD$Dscan/2, xmx=max(sD$Div.table$lon)+sD$Dscan/2, 
                                    ymn=min(sD$Div.table$lat)-sD$Dscan/2, ymx=max(sD$Div.table$lat)+sD$Dscan/2,resolution=sD$Dscan)
    ncell(r.div.refD)
    
    #plot div for ref Dist
    r.div.refD<-setValues(r.div.refD,sD$Div.table.cleaned$div.refD)
    {plot(r.div.refD,addfun=points(list_datageo[[i]]$lon,list_datageo[[i]]$lat,cex=0.2))
      title(paste(names(list_gen)[i],"div.refD",sep = " "),sub=paste0("Radius=",sD$Radius,"km, refD=",sD$refD,"km, ScanResol=",sD$Dscan,"?, Ncells=",sum(!is.na(sD$Div.table.cleaned$div.refD))))}
    
    #plot mean divergence
    r.div.mean<-setValues(r.div.mean,sD$Div.table$div.mean)
    {plot(r.div.mean,addfun=points(list_datageo[[i]]$lon,list_datageo[[i]]$lat,cex=0.2))
      title(paste(names(list_gen)[i],"div.mean",sep = " "),sub=paste0("Radius=",sD$Radius,"km, refD=",sD$refD,"km, ScanResol=",sD$Dscan,"?, Ncells=",sum(!is.na(sD$Div.table$div.mean))))}
    
    #
    plot(sD$Div.table.cleaned$subsample.size,sD$Div.table.cleaned$div.refD,ylim=c(0,1),title(paste(names(list_gen)[i]),sub=paste("N=", sum(!is.na(sD$Div.table.cleaned$div.refD))) ) ) #intrapolated divergence
    
    #Compare mean-div with interpolated div
    plot(sD$Div.table.cleaned$div.mean,sD$Div.table.cleaned$div.refD,ylim=c(0,1),title(paste(names(list_gen)[i]),sub=paste("N=", sum(!is.na(sD$Div.table.cleaned$div.refD))) ) ) #intrapolated divergence
    
    ##
    ## storing data across different parameters
    list_gen_dist[[paste("R",Radius,names(list_datageo)[i],sep='_')]]<-mHij
    list_r.div.refD[[paste("R",Radius,names(list_datageo)[i],sep='_')]]<-r.div.refD
    list_r.div.mean[[paste("R",Radius,names(list_datageo)[i],sep='_')]]<-r.div.mean
    list_sD[[paste("R",Radius,names(list_datageo)[i],sep='_')]]<-sD
    
  }
  
  dev.off()
  
}

str(list_sD)


list_geo_fut<-list_geo
list_datageo_fut<-list_datageo
list_gen_fut<-list_gen
list_gen_dist_fut<-list_gen_dist

##############
# SAVE IMAGE #
##############

save.image(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/futureDiv.Rdata")


