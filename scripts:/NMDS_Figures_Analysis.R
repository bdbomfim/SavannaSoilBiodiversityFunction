#vegetation analysis final
#ms PEA chapter2

library(vegan)
library(lattice)
library(permute)
library(MASS)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(patchwork)
citation()

#NMDS CAP 2 APPLIED SOIL ECOLOGY

#Upload seasonal fauna data
meso_cap2<-read.csv(file.choose())#20210502_FAUNA_NMDS_ASE
str(meso_cap2)
names(meso_cap2)

#data for analysis and figures
meso_cap2.1 = as.data.frame(sapply(meso_cap2[3:27], as.numeric))
str(meso_cap2.1)
#Plot names
rownames_cap2 <- meso_cap2$Season
rownames_cap2
groups<-unique(meso_cap2[0,3:27])
groups

#Environmental matrix seasonal
env_ss<-read.csv(file.choose())#20210503_Cap2_Jonas_Correlation
str(env_ss)
names(env_ss)
names(env_ss)[names(env_ss) == "Total_abundance"] <- "Abundance"
names(env_ss)[names(env_ss) == "Pielou"] <- "Evenness"
names(env_ss)[names(env_ss) == "Shannon"] <- "Diversity"
names(env_ss)[names(env_ss) == "Total_litterfall_g_m2_month"] <- "Total litterfall"
names(env_ss)[names(env_ss) == "Leaf_fall_g_m2_month"] <- "Leaf fall"
names(env_ss)[names(env_ss) == "Wood_fall_g_m2_month"] <- "Wood fall"
names(env_ss)[names(env_ss) == "Mean_monthly_rainfall_mm"] <- "Monthly rainfall"
names(env_ss)[names(env_ss) == "Leaf_k_at_720"] <- "Leaf_k"
names(env_ss)[names(env_ss) == "Wood_k_at_720"] <- "Wood_k"

#Seasonal environmental data frame
meso_cap2_env_b <- cbind(data.frame(env_ss[,c(3,5:8,11:13)]))
str(meso_cap2_env_b)

meso_cap2_env <- cbind(data.frame(meso_cap2[1:2],Plot_raw=Plot,env_ss[,c(3,5:8,11:13)]))
str(meso_cap2_env)

#Annual environmental data frame
env<-read.csv(file.choose())#20210503_Cap2_Env_Matrix
str(env)

#updating column names
names(env)[names(env) == "t"] <- "CECe"
names(env)[names(env) == "T"] <- "CECt"
names(env)[names(env) == "Mean_Richness"] <- "Richness"
names(env)[names(env) == "Annual_Abundance"] <- "Abundance"
names(env)[names(env) == "Mean_Diversity"] <- "Diversity"
names(env)[names(env) == "Mean_Evenness"] <- "Evenness"
names(env)[names(env) == "P.Rem"] <- "Prem"
names(env)[names(env) == "P"] <- "Pav"
names(env)[names(env) == "C_stock"] <- "SOCstock"
names(env)[names(env) == "H.Al"] <- "H+Al"
names(env)[names(env) == "OM"] <- "SOM"

#checking if names are alright
names(finalpca)
names(env)

#subsetting to include same final variables in PCA
env.1 = as.data.frame(sapply(env[,c(2:5,10:12,16,19:21,23)], as.numeric))
names(env.1)

##Load Annual Abundances matrix
yr_abund_matrix<-read.csv(file.choose())#20210503_FAUNA_NMDS_ANNUAL_ABUNDANCES
str(yr_abund_matrix)
names(yr_abund_matrix)
meso_cap2.1_an = as.data.frame(sapply(yr_abund_matrix[2:26], as.numeric))
str(meso_cap2.1_an)

# Community ecology variables####
H<-diversity(meso_cap2.1, index = "shannon")
H
S <- specnumber(meso_cap2.1) ## rowSums(BCI > 0) does the same...
S
J <- H/log(S)#Pielou's evenness J = H/log(S)
J
Plot<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
Plot
#Data frame for plotting
comm_data<-cbind(data.frame(meso_cap2$Season,meso_cap2$Plot,Plot,H,S,J))
comm_data

####Normalizing and running decostand####
head(meso_cap2.1)
new_mds <- vegdist(decostand(meso_cap2.1,"norm"), "bray")#distance matrix using bray curtis dissimilarity
new_mds

meso.mds_up2 <- metaMDS(new_mds,trymax=1000,k=2,autotransform = F)#species level based on dissimilarities
meso.mds_up2#Stress 0.06

p_env<-envfit(meso.mds_up2, meso_cap2_env_b,permutations=999)# this fits environmental vectors
p_env

p_spp.fit <- envfit(meso.mds_up2, meso_cap2.1, permutations = 999) # this fits species vectors
p_spp.fit

#PERMANOVA
adonis2(new_mds ~ meso_cap2$Season,permutations = 999,method="bray",by="margin")

#Preparing data for figures
site.scrs <- as.data.frame(scores(meso.mds_up2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scrs <- cbind(site.scrs, Season = meso_cap2_env$Season)
site.scrs <- cbind(site.scrs, Plot = meso_cap2_env$Plot_raw)
head(site.scrs)

species.scores <- as.data.frame(scores(p_spp.fit, display = "vectors"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores <- cbind(species.scores, Species = rownames(species.scores))
species.scores <- cbind(species.scores, pval = p_spp.fit$vectors$pvals)
sig.spp.scrs <- subset(species.scores, pval<=0.05)
head(species.scores)
head(sig.spp.scrs)
species.scores
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

env.scores <- as.data.frame(scores(p_env, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
env.scores <- cbind(env.scores, pval = p_env$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05
head(sig.env.scrs)
env.scores

nmds.plot.dune <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(Season), shape = factor(Season)), size = 4,alpha=0.8)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ scale_colour_manual(values = c("#fc4e51","#1dabe6"))+ theme(plot.title = element_text(size=18,hjust = 0.5))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Season", shape = "Season")+ # add legend labels for Management and Landuse
  theme(legend.position = "none", legend.text = element_text(size = 20), legend.title = element_text(size = 21), axis.text = element_text(size = 20,face="bold"),axis.title=element_text(size=22,face="bold")) # add legend at right of plot

nmds.plot.dune

#Adding soil epigeic fauna vectors
nmds.ssp<-nmds.plot.dune+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species),box.padding = 0.5,cex = 5.5, direction = "y", segment.size = 0.25,force=4,segment.alpha=0.6)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with soil fauna vectors")
nmds.ssp

#Adding seasonal env vectors
nmds.ssp.env<-nmds.ssp+
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lty=2,cex=0.8,colour = "#2F5BBA", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), box.padding = 2,cex = 5, direction = "x", segment.size = 0.4,colour = "#2F5BBA",force=4)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Soil epigeic fauna, litterfall and rainfall")
nmds.ssp.env

#Data frame for Figure 5a
MDS_xy<- as.data.frame(meso.mds_up2$points)
MDS_xy$Season<- meso_cap2$Season
Plot<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
Plot
MDS_xy$Plot<- as.factor(Plot)
meso_groups
str(MDS_xy)
head(meso_cap2)
MDS_xy2 <- rbind(data.frame(MDS_xy,Groups = meso_groups))
meso_groups
#Figure
p_nmds<-ggplot(MDS_xy, aes(MDS1, MDS2, shape = Season, colour = Season)) + geom_point(size=5,alpha=0.8) + theme_classic()+
  scale_colour_manual(values = c("#fc4e51","#1dabe6"))+
  theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),
  axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.title = element_text (size = 21), 
        legend.text = element_text (size = 20),legend.position=c(0.86,0.14),
        legend.background = element_rect(fill=alpha('transparent', 0.4)),legend.key=element_rect(fill=alpha('transparent', 0.4)))+labs(x = "NMDS1", y = "NMDS2")+ 
  annotate(geom="text", x=0.45, y=0.7, label="Stress = 0.06",color="black", size = 5)+theme(plot.title = element_text(size=18,hjust = 0.5))+
  stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  #annotate(geom="text", x=-0.7, y=0.7, label="a",color="black", size = 8,fontface="bold")+ 
  annotate(geom="text", x=0.38, y=0.62, label="PERMANOVA p=0.001",color="black", size = 5)+stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  labs(title = "Seasonal clustering")#+theme(axis.line = element_line(colour = "black", 

p_nmds



#Panel with both NMDS new Figure 5
p_nmds_all<-p_nmds+nmds.ssp.env+ plot_layout(ncol=2,widths=c(1,1))+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 22),plot.tag.position =c(0.05,1))
p_nmds_all

##Annual values####

new_mds_an <- vegdist(decostand(meso_cap2.1_an,"norm"), "bray")#distance matrix using bray curtis dissimilarity
new_mds_an

meso.mds_up2_an <- metaMDS(new_mds_an,trymax=1000,autotransform = F)#species level based on dissimilarities
meso.mds_up2_an#Stress 0.06

p_env_an<-envfit(meso.mds_up2_an, env.1,permutations=999)# this fits environmental vectors
p_env_an

p_spp.fit_an <- envfit(meso.mds_up2_an, meso_cap2.1_an, permutations = 999) # this fits species vectors
p_spp.fit_an

site.scrs_an <- as.data.frame(scores(meso.mds_up2_an))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#site.scrs_an <- cbind(site.scrs_an, Season = meso_cap2_env$Season)
site.scrs_an <- cbind(site.scrs_an, Plot = yr_abund_matrix$Plot)
site.scrs_an

species.scores_an <- as.data.frame(scores(p_spp.fit_an, display = "vectors"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores_an <- cbind(species.scores_an, Species = rownames(species.scores_an))
species.scores_an <- cbind(species.scores_an, pval = p_spp.fit_an$vectors$pvals)
sig.spp.scrs_an <- subset(species.scores_an, pval<0.05)
head(species.scores_an)
head(sig.spp.scrs_an)
species.scores_an
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

env.scores_an <- as.data.frame(scores(p_env_an, display = "vectors")) #extracts relevant scores from envifit
env.scores_an <- cbind(env.scores_an, env.variables_an = rownames(env.scores_an)) #and then gives them their names
env.scores_an <- cbind(env.scores_an, pval = p_env_an$vectors$pvals) # add pvalues to dataframe
env.scores_an
sig.env.scrs_an <- subset(env.scores_an, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores_an)
head(sig.env.scrs_an)

nmds.plot.dune_an <- ggplot(site.scrs_an, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2), shape=21,size = 2,alpha=0.8, color="black")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ theme(plot.title = element_text(size=18,hjust = 0.5))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  #labs(colour = "Season", shape = "Season")+ # add legend labels for Management and Landuse
  annotate(geom="text", x=0.6, y=0.75, label="Stress < 0.01",color="black", size = 5)+
  theme(legend.position = "none", legend.text = element_text(size = 20), 
        legend.title = element_text(size = 21), axis.text = element_text(size = 20,face="bold"),
        axis.title=element_text(size=22,face="bold")) # add legend at right of plot

nmds.plot.dune_an
sig.spp.scrs_an

#Adding species
nmds.ssp_an<-nmds.plot.dune_an+
  geom_segment(data = sig.spp.scrs_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3)+#geom_point(color = ifelse(dat2$car == "", "grey50", "red")) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs_an, aes(x=NMDS1, y=NMDS2, label = Species), box.padding = 0.8,cex = 5.5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundances")
nmds.ssp_an

#Adding all species
nmds.ssp_an_2<-nmds.plot.dune_an+
  geom_segment(data = species.scores_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2,colour=ifelse(pval<0.05, 'red', 'black')), arrow = arrow(length = unit(0.25, "cm")), colour = "grey50", lwd=0.3)+#geom_point(color = ifelse(dat2$car == "", "grey50", "red")) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = species.scores_an, aes(x=NMDS1, y=NMDS2, label = Species,colour=ifelse(pval<0.05, 'red', 'black')), box.padding = 0.8,cex = 5.5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundances")
nmds.ssp_an_2

#Adding env vectors
nmds.ssp.env_an<-nmds.ssp_an+
  geom_segment(data = sig.env.scrs_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lty=2,colour = "#2F5BBA", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.env.scrs_an, aes(x=NMDS1, y=NMDS2, label = env.variables_an), colour = "#2F5BBA",box.padding = 2,cex = 5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundance + decomposition")
nmds.ssp.env_an

#Adding all vectors that are in the PCA
nmds.ssp.env_an_2<-nmds.ssp_an+
  geom_segment(data = env.scores_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2,colour=ifelse(pval<0.05, 'darkred', 'grey10')), arrow = arrow(length = unit(0.25, "cm")), lty=2,lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = env.scores_an, aes(x=NMDS1, y=NMDS2, label = env.variables_an,colour=ifelse(pval<0.05, 'darkred', 'grey10')), box.padding = 0.6,cex = 5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundance + decomposition")
nmds.ssp.env_an_2


#Panel with both NMDS new Figure 5 HERE add all three - change lines to dashed when not significant
p_nmds_all_2<-p_nmds+nmds.ssp.env+pca_final+nmds.ssp.env_an+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 24,face="bold"),plot.tag.position =c(0.05,1))
p_nmds_all_2

#saving
ggsave(filename = "Fig5_Final_Panel_v7.png",
       plot = p_nmds_all_2, width = 20, height = 18, units = 'cm',
       scale = 2, dpi = 1000)

##PERMANOVA##

adonis2 (meso.1  ~ Littertrial+Net_tree_growth_m2.ha+wood.decay+Shannon + P.av.20.40 +P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="margin",data=env.4)


#END##
