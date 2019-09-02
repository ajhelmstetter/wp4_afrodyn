############################################################################
#Function to estimate diversity in a circular window of given Radius (in km)
#centered on latRef, lonRef (in degrees decimal) and accounting for 
#isolation by distance by estimating diversity as the expected divergence 
#between individuals at a given reference distance (refD, default refD = Radius/2).
#"mHij" is the matrix containing divergence values between all pairs of samples. 
#"datageo" must have columns named "lat" and "lon" with latitude and longitudes 
#of samples in degrees decimal, in the same order as for the "mHij" matrix.
#"geodist" is the distance object between individuals of "datageo". If it is
#not given when calling the function, it will be computed from "datageo" 
#but this is slower if the function is called repeatedly.
localDiv<-function(latRef,lonRef,Radius,refD=Radius/2,datageo,geodist=NULL,mHij){
  if(is.null(geodist)){
    print("Compute geodist")
    geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat))
  }
  
  #  dist_to_Ref=sqrt((datageo$lat-lat)^2+(datageo$lon-lon)^2)*40000/360 #distance in km btw Ref coordinate and each sample
  #earth great-circle distances of each sample to the ref coordinates
  dist_to_Ref <- 6371*acos(sin(latRef*pi/180)*sin(datageo$lat*pi/180) 
                           + cos(latRef*pi/180)*cos(datageo$lat*pi/180) 
                           * cos(datageo$lon*pi/180-lonRef*pi/180)) 
  
  keptind=dist_to_Ref<Radius 
  subsample=sum(keptind) #number of ind within circular window
  if(subsample>=3){
    subHij=as.dist(mHij[keptind,keptind])
    subgeodist=as.dist(as.matrix(geodist)[keptind,keptind])
    if(sum(subgeodist>0)>1){
      reg<-lm(subHij[subgeodist>0]~log(subgeodist[subgeodist>0]))
      HrefD<-reg$coefficients[1]+log(refD)*reg$coefficients[2]
    } else HrefD=NA
    Hmean<-mean(subHij)
    #  plot(x=log(subgeodist),y=subHij)
    #  points(y=HrefD,x=log(refD),col="red",cex=2)
    NshortD=sum(as.dist(subgeodist)<=refD)
    NhighD=sum(as.dist(subgeodist)>refD)
  }
  else HrefD=Hmean=NshortD=NhighD=NA
  if(subsample==2) Hmean=mHij[keptind,keptind][1,2]
  output=c(latRef,lonRef,HrefD,Hmean,subsample,NshortD,NhighD,Radius,refD)
  names(output)=c("lat","lon","div.refD","div.mean","subsample.size","Npairs.lowerD","Npairs.higherD","Radius(km)","Ref.Dist(km)")
  return(output)
}

compGen<-function(dataind,refDglobal){
#compute genetic distances between individuals for SNP or SSR data (unordered alleles)
Nind=nrow(dataind@tab) #number of individuals
Nloci=length(locNames(dataind)) #number of loci
Hijl<-array(dim=c(Nind,Nind,Nloci)) #array to store inter-individual genetic distances per locus (proportions of non identical pairs of alleles)
for(loc in 1:Nloci){ #loop over loci
  tab<-dataind[loc=loc]
  Hij=1-crossprod(t(tab@tab))/4  #proportion of non-identical pairs of alleles between diploid individuals
  Hij[is.na(Hij)]=mean(as.dist(Hij),na.rm=T) #replace missing values by overall means (between individuals only)
  Hijl[,,loc]=Hij

}
mHij=apply(Hijl,FUN="mean",MARGIN=c(1,2)) #pairwise genetic distances between individuals
return(mHij)
plot(x=log(geodist),y=as.dist(mHij))
#Regress mHij on log(geodist) to estimate divergence at a given ref distance
reg<-lm(as.dist(mHij)[geodist>0]~log(geodist[geodist>0]))
Div.global.refD<-reg$coefficients[1]+log(refDglobal)*reg$coefficients[2]
Div.global.refD
}

##
#Dscan=0.5
#Radius=100
#refD=10
#
#Scan local diversity around each node of a lon/lat grid of resolution Dscan in degrees
#NB: lat/lon scan order follows the order implicit in raster object (N->S and W->E)
scanDiv<-function(Dscan,Radius,refD,datageo,geodist){
Div.table=NULL
for(lat in Dscan*(ceiling(max(datageo$lat)/Dscan):floor(min(datageo$lat)/Dscan))){
  for(lon in Dscan*(floor(min(datageo$lon)/Dscan):ceiling(max(datageo$lon)/Dscan))){
    print(c(lat,lon))
    div=localDiv(latRef=lat,lonRef=lon,Radius=Radius,refD=refD,datageo=datageo,geodist=geodist,mHij=mHij)
    if(is.null(Div.table)) Div.table=div else Div.table<-rbind(Div.table,div)
  }
}
row.names(Div.table)<- 1:length(Div.table[,1])
Div.table=as.data.frame(Div.table)
#clean the table to keep only estimates of divergence at ref distance based on interpolation and not extrapolation)
Div.table.cleaned=Div.table
Div.table.cleaned$div.refD[Div.table$Npairs.lowerD==0 | Div.table$Npairs.higherD==0]=NA
#short table limited to cells with estimates of div.mean
Div.table.short=Div.table[!is.na(Div.table$div.mean),]

newList <- list("Div.table" = Div.table,
                "Div.table.cleaned" = Div.table.cleaned,
                "Dscan" = Dscan,
                "Radius" = Radius,
                "refD" = refD)

return(newList)

#Check if estimates depends on sample size within local window
plot(Div.table$subsample.size,Div.table$div.mean,ylim=c(0,1),title(sub=paste("N=", sum(!is.na(Div.table$div.mean))) ) ) #mean divergence
plot(Div.table$subsample.size,Div.table$div.refD,ylim=c(0,1),title(sub=paste("N=", sum(!is.na(Div.table$div.refD))) ) ) #intra+extra-polated divergence
plot(Div.table.cleaned$subsample.size,Div.table.cleaned$div.refD,ylim=c(0,1),title(sub=paste("N=", sum(!is.na(Div.table.cleaned$div.refD))) ) ) #intrapolated divergence
#Compare mean-div with interpolated div
plot(Div.table.cleaned$div.mean,Div.table.cleaned$div.refD,ylim=c(0,1),title(sub=paste("N=", sum(!is.na(Div.table.cleaned$div.refD))) ) ) #intrapolated divergence
}


