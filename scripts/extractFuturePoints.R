setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")

load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/results/currentDiv.Rdata")
load(file="~/Dropbox/projects/AJH_WP4_AFRODYN/data/rasters.RData")
load("~/Dropbox/projects/AJH_WP4_AFRODYN/data/wp4.Rdata")

names(list_r.div.mean)
list_layers

library(raster)

#resample(list_layers[[1]]$loss_win_raster,list_r.div.mean$R_100_aa)

list_layers_mod<-list(list_layers$`Annickia affinis`,
                      list_layers$`Anonidium mannii`,
                      list_layers$`Friesodielsia enghiana`,
                      list_layers$`Greenwayodendron suaveolens`,
                      list_layers$`Podococcus barteri`,
                      list_layers$`Podococcus acaulis`,
                      list_layers$`Sclerosperma mannii`,
                      list_layers$`Distemonanthus benthamianus`,
                      list_layers$`Barteria fistulosa`,
                      list_layers$`Erythrophleum ivorense`,
                      list_layers$`Erythrophleum suaveolens`,
                      list_layers$`Aucoumea klaineana`,
                      list_layers$`Greenwayodendron suaveolens`,
                      list_layers$`Santiria trimera`,
                      list_layers$`Santiria trimera`,
                      list_layers$`Santiria trimera`,
                      list_layers$`Santiria trimera`,
                      list_layers$`Scorodophloeus zenkeri`,
                      list_layers$`Terminalia superba`
                      )
             
   
list_datageo_fut<-list()
list_gen_fut<-list()
list_geo_fut<-list()

#problem with 001/011 numbers
indNames(list_gen$bf)<-as.character(as.numeric(indNames(list_gen$bf)))
indNames(list_gen$ak)<-as.character(as.numeric(indNames(list_gen$ak)))

list_datageo$bf<-data.frame(list_datageo$bf$lon,list_datageo$bf$lat)
list_datageo$ak<-data.frame(list_datageo$ak$lon,list_datageo$ak$lat)
list_datageo$gs_ms<-data.frame(list_datageo$gs_ms$lon,list_datageo$gs_ms$lat)
list_datageo$s_ebo<-data.frame(list_datageo$s_ebo$lon,list_datageo$s_ebo$lat)
list_datageo$s_obov<-data.frame(list_datageo$s_obov$lon,list_datageo$s_obov$lat)
list_datageo$s_tr_ca<-data.frame(list_datageo$s_tr_ca$lon,list_datageo$s_tr_ca$lat)
list_datageo$s_tr_wa<-data.frame(list_datageo$s_obov$lon,list_datageo$s_tr_wa$lat)
list_datageo$s_zen<-data.frame(list_datageo$s_zen$lon,list_datageo$s_zen$lat)
list_datageo$t_sup<-data.frame(list_datageo$t_sup$lon,list_datageo$t_sup$lat)

colnames(list_datageo$bf)<-c("lon","lat")
colnames(list_datageo$ak)<-c("lon","lat")
colnames(list_datageo$gs_ms)<-c("lon","lat")
colnames(list_datageo$s_ebo)<-c("lon","lat")
colnames(list_datageo$s_obov)<-c("lon","lat")
colnames(list_datageo$s_tr_ca)<-c("lon","lat")
colnames(list_datageo$s_tr_wa)<-c("lon","lat")
colnames(list_datageo$s_zen)<-c("lon","lat")
colnames(list_datageo$t_sup)<-c("lon","lat")

rownames(list_datageo$gs_ms)<-indNames(list_gen$gs_ms)
rownames(list_datageo$s_ebo)<-indNames(list_gen$s_ebo)
rownames(list_datageo$s_obov)<-indNames(list_gen$s_obov)
rownames(list_datageo$s_tr_ca)<-indNames(list_gen$s_tr_ca)
rownames(list_datageo$s_tr_wa)<-indNames(list_gen$s_tr_wa)
rownames(list_datageo$s_zen)<-indNames(list_gen$s_zen)
rownames(list_datageo$t_sup)<-indNames(list_gen$t_sup)

for(i in 1:length(list_datageo)){
  #get individuals that are in non-declining areas 
  loss<-extract(list_layers_mod[[i]]$loss_win_raster,list_datageo[[i]])
  loss[is.na(loss)] <- 0
  tmp<-cbind(list_datageo[[i]],loss)
  tmp2<-tmp[tmp$loss!=-1,]
  list_datageo_fut[[i]]<-tmp2[,1:2]
  
  #subset genind objects
  list_gen_fut[[i]]<-list_gen[[i]][rownames(list_datageo_fut[[i]]), ]
  
  #get geographic distances
  list_geo_fut[[i]]<-as.dist(geod(lon=list_datageo_fut[[i]]$lon,lat=list_datageo_fut[[i]]$lat)) 
  
  
}

remove(list_layers)
remove(list_layers_mod)

save.image("data/datageo_fut.Rdata")
