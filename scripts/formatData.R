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

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Distemonanthus benthamianus/inK=5_250814(approx data Demenou 2017).txt",header = T,row.names=1)

###
#Geographic data
###

datageo<-data[,2:3]
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-3,-15)],sep="-",NA.char="0",type = "codom")

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


############################################################################################################
############################################################################################################
#####
# Erythrophleum ivorense
#####

data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/inSpag-Eivo.csv",header = T,row.names=1)

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
list_geo[['ei']]<-geodist
list_gen[['ei']]<-dataind
list_datageo[['ei']]<-datageo


############################################################################################################
############################################################################################################
#####
# Erythrophleum suaveolens (Ouest et centrale)
#####

# individual OH1579 had a duplicate, remove the one with missing data

data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/inSpag-EsuaAC-EsuaW.csv",header = T,row.names=1)

#remove genetic pool column
data<-data[,-1]

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
list_geo[['es']]<-geodist
list_gen[['es']]<-dataind
list_datageo[['es']]<-datageo


############################################################################################################
############################################################################################################
#####
# Aukoumea klaineana
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Aukoumea klaineana/inAukoumea542ind(data_Born2011MolEcol).txt",header = T)

###
#Geographic data
###

datageo<-data[,c(4,3)]
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-4,-19)],ncode=3,NA.char="0",type = "codom")

#append to list
list_geo[['ak']]<-geodist
list_gen[['ak']]<-dataind
list_datageo[['ak']]<-datageo

############################################################################################################
############################################################################################################
#####
# Greenwayodendron suaveolens MICROSAT
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Greenwayodendron/inSpag9popGeo.txt",header = T,row.names=1)

###
#Geographic data
###

datageo<-data[,3:2]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-3,-12,-13)],ncode=3,NA.char="0",type = "codom")

#append to list
list_geo[['gs_ms']]<-geodist
list_gen[['gs_ms']]<-dataind
list_datageo[['gs_ms']]<-datageo

############################################################################################################
############################################################################################################
#####
# Santiria ebo
#####

data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Santiria/s_ebo.csv",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,2:1]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-2)],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['s_ebo']]<-geodist
list_gen[['s_ebo']]<-dataind
list_datageo[['s_ebo']]<-datageo


############################################################################################################
############################################################################################################
#####
# Santiria obovata
#####

#GiD2255 duplicated so 1 removed


data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Santiria/s_obovata.csv",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,2:1]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-2)],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['s_obov']]<-geodist
list_gen[['s_obov']]<-dataind
list_datageo[['s_obov']]<-datageo


############################################################################################################
############################################################################################################
#####
# Santiria trimera central Africa
#####

data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Santiria/s_trim_CA.csv",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,2:1]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-2)],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['s_tr_ca']]<-geodist
list_gen[['s_tr_ca']]<-dataind
list_datageo[['s_tr_ca']]<-datageo

############################################################################################################
############################################################################################################
#####
# Santiria trimera west Africa
#####

data<-read.csv("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Santiria/s_trim_WA.csv",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,2:1]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-2)],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['s_tr_wa']]<-geodist
list_gen[['s_tr_wa']]<-dataind
list_datageo[['s_tr_wa']]<-datageo


############################################################################################################
############################################################################################################
#####
# Scorodophleus zenkeri
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Scorodophleus zenkeri/inScoro466ind6clusters.txt",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,3:2]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-3)],ncode=3,NA.char="0",type = "codom")

#append to list
list_geo[['s_zen']]<-geodist
list_gen[['s_zen']]<-dataind
list_datageo[['s_zen']]<-datageo

############################################################################################################
############################################################################################################
#####
# Terminalia superba
#####

data<-read.table("~/Dropbox/projects/AJH_WP4_AFRODYN/data/Terminalia superba/Inspagedi271-5pools.txt",header = T,row.names=1)

###
#Geographic data
###
str(data)
datageo<-data[,3:2]
datageo
geodist<-as.dist(geod(lon=datageo$lon,lat=datageo$lat)) #great-circle distances (km)

###
# Microsats
###

#extract genotypic data in genind format (adegenet). sep defines the character separating alleles within locus
dataind<-df2genind(data[,c(-1:-3,-18)],sep="-",NA.char="0",type = "codom")

#append to list
list_geo[['t_sup']]<-geodist
list_gen[['t_sup']]<-dataind
list_datageo[['t_sup']]<-datageo

##############
# SAVE IMAGE #
##############

save.image(file="~/Dropbox/projects/AJH_WP4_AFRODYN/data/wp4.Rdata")
