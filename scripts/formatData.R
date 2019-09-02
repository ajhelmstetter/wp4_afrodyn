rm(list=ls())
library(adegenet)
library(pegas)
library(vcfR)
#library(PopGenome)
#library(poppr)

setwd("~/Dropbox/projects/AJH_WP4_AFRODYN/")
source("scripts/localDiv.R")

############################################################################################################
############################################################################################################
#####
# Annickia affinis
#####

###### Geographic data

id_ind<-read.csv("~/Dropbox/projects/AJH_AFRODYN/annickia/ID_index.csv")
id_vou<-read.csv("~/Dropbox/projects/AJH_AFRODYN/annickia/ID_voucher.csv")

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

###### SNP

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/annickia/anni_filtered.sub2.vcf",checkFile = T, convertNA = T) #read in all data

# convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order

datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]

data.frame(indNames(dataind),rownames(datageo))

#create list
list_geo<-list('aa'=geodist)
list_gen<-list('aa'=dataind)
list_datageo<-list('aa'=datageo)

############################################################################################################
############################################################################################################
#####
# Anonidium mannii
#####

###
#Geographic data
###

id_ind<-read.csv("~/Dropbox/projects/AJH_AFRODYN/anonidium/ID_index.csv")
id_vou<-read.csv("~/Dropbox/projects/AJH_AFRODYN/anonidium/ID_voucher.csv")
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
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/anonidium/anon_filtered.sub2.vcf",checkFile = T, convertNA = T) #read in all data

### convert to genlight
dataind <- vcfR2genind(vcf)

indNames(dataind)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]

data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['am']]<-geodist
list_gen[['am']]<-dataind
list_datageo[['am']]<-datageo

############################################################################################################
############################################################################################################
#####
# Monanthotaxis enghiana
#####

###
#Geographic data
###

id_ind<-read.csv("~/Dropbox/projects/AJH_AFRODYN/fries/ID_index.csv")
id_vou<-read.csv("~/Dropbox/projects/AJH_AFRODYN/fries/ID_voucher.csv")

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])

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
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/fries/fries_filtered.sub2.vcf") #read in all data

### convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]
data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['me']]<-geodist
list_gen[['me']]<-dataind
list_datageo[['me']]<-datageo

############################################################################################################
############################################################################################################
#####
# Greenwayodendron suaveolens
#####

###
#Geographic data
###

id_ind<-read.csv("~/Dropbox/projects/AJH_AFRODYN/greenwayodendron/ID_index.csv")
id_vou<-read.csv("~/Dropbox/projects/AJH_AFRODYN/greenwayodendron/ID_voucher.csv")

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
foo<-foo[order(foo$index),]
foo<-foo[c(1:34,54:length(foo$index)),]

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
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/greenwayodendron/green_filtered.sub2.vcf") #read in all data

### convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]
data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['gs']]<-geodist
list_gen[['gs']]<-dataind
list_datageo[['gs']]<-datageo

############################################################################################################
############################################################################################################
#####
# Podococcus barteri
#####

###
#Geographic data
###

foo<-read.csv("~/Dropbox/projects/AJH_AFRODYN/podo/foo.csv")

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
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/podo/podo_b_filtered.sub2.vcf") #read in all data

### convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]
data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['pb']]<-geodist
list_gen[['pb']]<-dataind
list_datageo[['pb']]<-datageo

############################################################################################################
############################################################################################################
#####
# Podococcus acaulis
#####

###
#Geographic data
###
foo<-read.csv("~/Dropbox/projects/AJH_AFRODYN/podo_a/foo.csv")

#extract geographic data
#name | lon | lat
datageo<-data.frame(as.numeric(foo$ddlong),as.numeric(foo$ddlat))
rownames(datageo)<-foo$index
datageo
colnames(datageo)<-c("lon","lat")
#calculate genetic distances
#great-circle distances (km)
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) 

###
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/podo_a/podo_a_filtered.sub2.vcf") #read in all data

### convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]
data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['pa']]<-geodist
list_gen[['pa']]<-dataind
list_datageo[['pa']]<-datageo

############################################################################################################
############################################################################################################
#####
# Sclerosperma mannii
#####

###
#Geographic data
###

id_ind<-read.csv("~/Dropbox/projects/AJH_AFRODYN/sclero/ID_index.csv")
id_vou<-read.csv("~/Dropbox/projects/AJH_AFRODYN/sclero/ID_voucher.csv")

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])

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
# SNP
###

vcf <- read.vcfR("~/Dropbox/projects/AJH_AFRODYN/sclero/sclero_filtered.sub2.vcf") #read in all data

### convert to genind
dataind <- vcfR2genind(vcf)

## Check geo and genetic data are in same order
datageo<-datageo[order(match(rownames(datageo),indNames(dataind))),]
data.frame(indNames(dataind),rownames(datageo))

#append to list
list_geo[['sm']]<-geodist
list_gen[['sm']]<-dataind
list_datageo[['sm']]<-datageo

############################################################################################################
############################################################################################################
#####
# Distemonanthus	benthaniaus
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Distem381ind.txt",header = T,row.names=1)

###
#Geographic data
###

datageo<-data[,1:2]
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,-1:-2],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['db']]<-geodist
list_gen[['db']]<-dataind
list_datageo[['db']]<-datageo


############################################################################################################
############################################################################################################
#####
# Barteria fistulosa
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Appendix1.txt",header = T)

###
#Geographic data
###

datageo<-data[,2:3]
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,-1:-3],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['bf']]<-geodist
list_gen[['bf']]<-dataind
list_datageo[['bf']]<-datageo




##############
# SAVE IMAGE #
##############

save.image(file="~/Dropbox/projects/AJH_WP4_AFRODYN/data/wp4.Rdata")
