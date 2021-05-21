##MONTE CARLO TEST AND ANNUAL FAUNA DATA

library(sp)
library(spdep)
library(tmap)
library(ade4)
library(geosphere)
library(vegan)
load(url("http://github.com/mgimond/Spatial/raw/master/Data/moransI.RData"))
class(s1)
xy <- sp::coordinates(s1)
xy

#Upload annual fauna data
ann_meso_cap2<- read.csv(file.choose())#20210503_FAUNA_NMDS_ANNUAL_ABUNDANCES

#Examine data
names(ann_meso_cap2)
ann_meso_cap2<-as.data.frame(sapply(ann_meso_cap2[1:31], as.numeric))
str(ann_meso_cap2)
ann_meso_cap2$Plot=as.factor(ann_meso_cap2$Plot)

#Mantel test for spatial autocorrelation
geo<-cbind(ann_meso_cap2$Longitude,ann_meso_cap2$Latitude)
plot.dists <- dist(geo)
#geographic data frame - haversine distance 
plot.dists <- distm(geo, fun = distHaversine)
dist.geo = as.dist(plot.dists)

fauna.dists <- dist(cbind(ann_meso_cap2$Total_Abundance,ann_meso_cap2$Mean_Diversity),method="euclidean")
fauna.dists2 <- vegdist(ann_meso_cap2[2:27],method="bray")

as.matrix(plot.dists)
as.matrix(fauna.dists)
as.matrix(fauna.dists2)

#abundance vs geographic
mantel.rtest(plot.dists, fauna.dists, nrepet = 999)#r=0.09

#abundance vs geographic
mantel.rtest(plot.dists, fauna.dists2, nrepet = 999)#r=0.03

##END###
