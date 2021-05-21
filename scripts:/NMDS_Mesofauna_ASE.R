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
names(env_ss)[names(env_ss) == "Shannon"] <- "Diversity"
names(env_ss)[names(env_ss) == "Total_litterfall_g_m2_month"] <- "Tot_Litterfall"
names(env_ss)[names(env_ss) == "Leaf_fall_g_m2_month"] <- "Leaf_fall"
names(env_ss)[names(env_ss) == "Wood_fall_g_m2_month"] <- "Wood_fall"
names(env_ss)[names(env_ss) == "Mean_monthly_rainfall_mm"] <- "Monthly_rainfall"

#by season environmental data frame
meso_cap2_env <- cbind(data.frame(meso_cap2[1:2],Plot_raw=Plot,env_ss[,c(6:8,11:13)]))
str(meso_cap2_env)

#Environmental matrix annual
env<-read.csv(file.choose())#20210503_Cap2_Env_Matrix
str(env)
names(env)
env.1 = as.data.frame(sapply(env[,c(2:13,15:19)], as.numeric))
str(env.1)
env.2=as.data.frame(sapply(env[,c(2:10,15:16,18:19)], as.numeric))
env.2
env.3=env[,c(1:2,5,8,10,14,17,24,31,44:45)]
env.3
env.4=env[,c(1:2,5,8,10,14,17,31,44:45)]
env.4
env.5=env[,c(1:2,5,10,14,17,31,44:45)]
env.5

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

p_env<-envfit(meso.mds_up2, meso_cap2_env,permutations=999)# this fits environmental vectors
p_env

p_spp.fit <- envfit(meso.mds_up2, meso_cap2.1, permutations = 999) # this fits species vectors
p_spp.fit

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
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

env.scores <- as.data.frame(scores(p_env, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
env.scores <- cbind(env.scores, pval = p_env$vectors$pvals) # add pvalues to dataframe
sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05
head(sig.env.scrs)

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

#Data frame
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
  annotate(geom="text", x=0.45, y=0.7, label="Stress = 0.06",color="black", size = 6)+theme(plot.title = element_text(size=18,hjust = 0.5))+
  stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  #annotate(geom="text", x=-0.7, y=0.7, label="a",color="black", size = 8,fontface="bold")+ 
  annotate(geom="text", x=0.3, y=0.62, label="PERMANOVA p=0.001",color="black", size = 6)+stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  labs(title = "Seasonal clustering")#+theme(axis.line = element_line(colour = "black", 

p_nmds

#Panel with both NMDS new Figure 5
p_nmds_all<-p_nmds+nmds.ssp.env+ plot_layout(ncol=2,widths=c(1,1))+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 22),plot.tag.position =c(0.05,1))
p_nmds_all

#Saving figure in high res new Figure 5
ggsave(filename = "Fig5_Final_NMDS_v7.png",
       plot = p_nmds_all, width = 18, height = 9, units = 'cm',
       scale = 2, dpi = 1000)

##PREPARE THIRD PANEL WITH SOIL VECTORS - check coinertia analysis####

#1 - run new NMDS WITH annual abundances as one matrix and soil and litter variables as another matrix
#dataframes are meso_anab and env
#2 - run envfit and adonis2
#3 - fit vectors for faunal groups and decomposition and soil variables

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
#species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

env.scores_an <- as.data.frame(scores(p_env_an, display = "vectors")) #extracts relevant scores from envifit
env.scores_an <- cbind(env.scores_an, env.variables_an = rownames(env.scores_an)) #and then gives them their names
env.scores_an <- cbind(env.scores_an, pval = p_env_an$vectors$pvals) # add pvalues to dataframe
env.scores_an
sig.env.scrs_an <- subset(env.scores_an, pval<=0.05) #subset data to show variables significant at 0.05
head(env.scores_an)
head(sig.env.scrs_an)

nmds.plot.dune_an <- ggplot(site.scrs_an, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2), shape=21,size = 3,alpha=0.8, color="black")+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ theme(plot.title = element_text(size=18,hjust = 0.5))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  #labs(colour = "Season", shape = "Season")+ # add legend labels for Management and Landuse
  theme(legend.position = "none", legend.text = element_text(size = 20), 
        legend.title = element_text(size = 21), axis.text = element_text(size = 20,face="bold"),
        axis.title=element_text(size=22,face="bold")) # add legend at right of plot

nmds.plot.dune_an
sig.spp.scrs_an

#Adding species
nmds.ssp_an<-nmds.plot.dune_an+
  geom_segment(data = sig.spp.scrs_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey50", lwd=0.3)+#geom_point(color = ifelse(dat2$car == "", "grey50", "red")) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs_an, aes(x=NMDS1, y=NMDS2, label = Species), box.padding = 0.8,cex = 5.5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundances")
nmds.ssp_an

#Adding env vectors
nmds.ssp.env_an<-nmds.ssp_an+
  geom_segment(data = sig.env.scrs_an, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lty=2,colour = "#2F5BBA", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.env.scrs_an, aes(x=NMDS1, y=NMDS2, label = env.variables_an), colour = "#2F5BBA",box.padding = 2,cex = 5, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Annual abundances and litter decomposition")
nmds.ssp.env_an

#Panel with both NMDS new Figure 5 HERE add all three - change lines to dashed when not significant
p_nmds_all_2<-p_nmds+nmds.ssp.env+nmds.ssp.env_an+plot_annotation(tag_levels = 'a')& 
  theme(plot.tag = element_text(size = 22),plot.tag.position =c(0.05,1))
p_nmds_all_2
(p1 | p2)/(p3 | p4)

#saving
ggsave(filename = "Fig5_NMDS_All_env_v8.png",
       plot = p_nmds_all_2, width = 24, height = 9, units = 'cm',
       scale = 2, dpi = 1000)

#Data frame
MDS_xy_an<- as.data.frame(meso.mds_up2$points)
MDS_xy$Season<- meso_cap2$Season
Plot<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
Plot
MDS_xy$Plot<- as.factor(Plot)
meso_groups
str(MDS_xy)
head(meso_cap2)
MDS_xy2 <- rbind(data.frame(MDS_xy,Groups = meso_groups))

#Figure
p_nmds<-ggplot(MDS_xy, aes(MDS1, MDS2, shape = Season, colour = Season)) + geom_point(size=5,alpha=0.8) + theme_classic()+
  scale_colour_manual(values = c("#fc4e51","#1dabe6"))+
  theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),
        axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.title = element_text (size = 21), 
        legend.text = element_text (size = 20))+labs(x = "NMDS1", y = "NMDS2")+ 
  annotate(geom="text", x=0.3, y=0.7, label="Stress = 0.06",color="black", size = 6)+theme(plot.title = element_text(size=18,hjust = 0.5))+
  stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  #annotate(geom="text", x=-0.7, y=0.7, label="a",color="black", size = 8,fontface="bold")+ 
  annotate(geom="text", x=0.3, y=0.65, label="PERMANOVA p=0.001",color="black", size = 6)+stat_ellipse(aes(colour=Season),type="t",linetype=2)+
  labs(title = "NMDS clustering")#+theme(axis.line = element_line(colour = "black", 

p_nmds


#### Envfit####

en_annual_v0<-envfit(meso_anab, env, permutations = 999, na.rm = TRUE)
en_annual_v0

en_annual<-envfit(meso_anab, env.1, permutations = 999, na.rm = TRUE)
en_annual
#K highly significant and P marginally significant

en_annual_2<-adonis2 (meso_anab  ~ P+K, permutations = 999,method="bray",by="margin",data=env.1)
en_annual_2

en_annual_v2<-envfit(meso_anab, env.1, permutations = 999, na.rm = TRUE)
en_annual

en = envfit(meso.mds_up2, env.1, permutations = 999, na.rm = TRUE)
en$factors

en2 <- envfit(meso.mds_up2, env.2, permutations = 999, na.rm = TRUE)
str(en2)
names(env.2)

en3 <- envfit(meso.mds_up2, env.3, permutations = 999, na.rm = TRUE)
en3
names(env.2)

en4 <- envfit(meso.mds_up2, env.4, permutations = 999, na.rm = TRUE)
en4

en5 <- envfit(meso.mds_up2, env.5, permutations = 999, na.rm = TRUE)
en5

plot(meso.mds_up2)
plot(en2)

en_coord_cont = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)
en_coord_cat = as.data.frame(scores(en2, "factors")) * ordiArrowMul(en2)

data.scores = as.data.frame(scores(meso.mds_up2))
data.scores$Stand <- meso$Stand
data.scores$LitterTrial <- meso$Litter.trial

gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LitterTrial, shape=Stand), size = 5)+ theme_replace()+scale_colour_manual(values = c("#fc4e51","#1dabe6"))+theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18))+labs(x = "NMDS1", y = "NMDS2")+ annotate(geom="text", x=-0.25, y=-0.22, label="Stress = 0.09",
                                                                                                                                  color="black", size = 6)

gg

gg2 <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LitterTrial, shape=Stand), size = 5)+ theme_replace()+scale_colour_manual(values = c("#fc4e51","#1dabe6"))+theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18))+labs(x = "NMDS1", y = "NMDS2")+ annotate(geom="text", x=-0.25, y=-0.22, label="Stress = 0.09",
                                                                                                                                  color="black", size = 6)+ 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont))

gg2

en_coord_cont3 = as.data.frame(scores(en3, "vectors")) * ordiArrowMul(en3)
en_coord_cat3 = as.data.frame(scores(en3, "factors")) * ordiArrowMul(en3)

en_coord_cont4 = as.data.frame(scores(en4, "vectors")) * ordiArrowMul(en4)
en_coord_cat4 = as.data.frame(scores(en4, "factors")) * ordiArrowMul(en4)

en_coord_cont5 = as.data.frame(scores(en5, "vectors")) * ordiArrowMul(en5)
en_coord_cat5 = as.data.frame(scores(en5, "factors")) * ordiArrowMul(en5)


data.scores = as.data.frame(scores(meso.mds_up2))
data.scores$Stand <- meso$Stand
data.scores$LitterTrial <- meso$Litter.trial

gg3 <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LitterTrial, shape=Stand), size = 5)+ theme_replace()+scale_colour_manual(values = c("#fc4e51","#1dabe6"))+theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18))+labs(x = "NMDS1", y = "NMDS2")+ annotate(geom="text", x=-0.25, y=-0.22, label="Stress = 0.09",
                                                                                                                                  color="black", size = 6)+ geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                                                                                                                                                                         data = en_coord_cont3, size =1, alpha = 0.5, colour = "grey30") +
  #geom_point(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2), 
             #shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2+0.04), 
            #label = row.names(en_coord_cat3), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont3, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont3), size=6,position = position_nudge(y = 0.01,x=0.02))

gg3+stat_ellipse(aes(colour=LitterTrial),type="t",linetype=2,size=1.4,alpha=0.6)

gg4 <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LitterTrial, shape=Stand), size = 5)+ theme_replace()+scale_colour_manual(values = c("#fc4e51","#1dabe6"))+theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18))+labs(x = "NMDS1", y = "NMDS2")+ annotate(geom="text", x=-0.25, y=-0.22, label="Stress = 0.09",
                                                                                                                                  color="black", size = 6)+ geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                                                                                                                                                                         data = en_coord_cont4, size =1, alpha = 0.5, colour = "grey30") +
  #geom_point(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2), 
  #shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2+0.04), 
  #label = row.names(en_coord_cat3), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont4, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont4), size=6,position = position_nudge(y = 0.01,x=0.01))

gg4+stat_ellipse(aes(colour=LitterTrial),type="t",linetype=2)

gg5 <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = LitterTrial, shape=Stand), size = 5)+ theme_replace()+scale_colour_manual(values = c("#fc4e51","#1dabe6"))+theme(axis.title=element_text(size=22,face="bold"),axis.text.y=element_text(size=20, face="bold"),axis.text.x = element_text(angle=0, hjust=0.5,size=20, face="bold"))+
  theme(legend.title = element_text (size = 20), legend.text = element_text (size = 18))+labs(x = "NMDS1", y = "NMDS2")+ annotate(geom="text", x=-0.25, y=-0.22, label="Stress = 0.09",
                                                                                                                                  color="black", size = 7)+ geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                                                                                                                                                                         data = en_coord_cont5, size =1, alpha = 0.5, colour = "grey30") +
  #geom_point(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2), 
  #shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  #geom_text(data = en_coord_cat3, aes(x = NMDS1, y = NMDS2+0.04), 
  #label = row.names(en_coord_cat3), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont5, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont5), size=6,position = position_nudge(y = 0.01,x=0.02))

gg5+stat_ellipse(aes(colour=LitterTrial),type="t",linetype=2,size=1.4)

#### RUN PERMANOVA ####

## Check assumptions as in here https://chrischizinski.github.io/rstats/adonis/

data(dune)
str(dune)
data(dune.env)
str(dune.env)

names(env.4)

adonis2 (meso.1  ~ Littertrial+Net_tree_growth_m2.ha+wood.decay+Shannon + P.av.20.40 +P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="margin",data=env.4)

adonis2 (meso.1  ~ Littertrial+Net_tree_growth_m2.ha+wood.decay+Shannon + P.av.20.40 +P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="terms",data=env.4)

adonis2 (meso.1  ~ Littertrial+Net_tree_growth_m2.ha+wood.decay+ P.av.20.40 +P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="terms",data=env.4)

adonis2 (meso.1  ~ Littertrial+wood.decay+Shannon +P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="margin",data=env.4)

adonis2 (meso.1  ~ Littertrial+wood.decay+P.rem.0.20 + Vol.moisture.0.20, permutations = 999,method="bray",by="margin",data=env.4)


#### Sample NMDS Plotting ####

## Sample 0

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Host = metadata$Host_Taxon, Location = metadata$Location)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Location)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "NMDS Plot")

## More here http://geoffreyzahn.com/nmds-example/

##Sample 1

# calculate distance for NMDS
sol <- metaMDS(dune)
# Create meta data for grouping
MyMeta = data.frame(
  sites = c(2,13,4,16,6,1,8,5,17,15,10,11,9,18,3,20,14,19,12,7),
  amt = c("hi", "hi", "hi", "md", "lo", "hi", "hi", "lo", "md", "md", "lo", 
          "lo", "hi", "lo", "hi", "md", "md", "lo", "hi", "lo"),
  row.names = "sites")
# plot NMDS using basic plot function and color points by "amt" from MyMeta
plot(sol$points, col = MyMeta$amt)
# draw dispersion ellipses around data points
ordiellipse(sol, MyMeta$amt, display = "sites", kind = "sd", label = T)

# same in ggplot2
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(data = MyMeta, color = MyMeta$amt))

ordiellipse(sol, MyMeta$amt, display = "sites", kind = "se", conf = 0.95, label = T)

##https://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo

##Sample 2

data.scores <- as.data.frame(scores(meso.mds_up2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site<-rownames(data.scores)   # create a column of site names, from the rownames of data.scores
data.scores$grp <- rownames  #  add the grp variable created earlier
head(data.scores)

species.scores <- as.data.frame(scores(meso.mds_up2, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)


## Sample 3

plot(dune.mds, type = "n") #displays empty ordination space
points(dune.mds, display = "sites", pch = c(16, 8, 17, 11) [as.numeric(dune.env$Management)], col = c("blue", "orange", "black") [as.numeric(dune.env$Use)]) # displays site points where symbols (pch) are different management options and colour (col) are different land uses
legend("topright", legend = c(levels(dune.env$Management), levels(dune.env$Use)), pch = c(16, 8, 17, 11, 16, 16, 16), col = c("black","black","black","black","blue", "orange", "black"), bty = "n", cex = 1) # displays symbol and colour legend
legend("topleft", "stress = 0.118", bty = "n", cex = 1) # displays legend text of stress value 

nmds.plot.dune <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$Management), shape = factor(site.scrs$Landuse)), size = 2)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  labs(colour = "Management", shape = "Landuse")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

nmds.plot.dune + labs(title = "Basic ordination plot") #displays plot

nmds.plot.dune+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")

#Plotting ENVIRONMENTAL VARIABLES

nmds.plot.dune+
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables), cex = 4, direction = "both", segment.size = 0.25)+ #add labels for env variables
  labs(title="Ordination with environmental vectors")

### MORE here https://www.rpubs.com/RGrieger/545184

##End Sample 3

sum(meso$FOR)
diversity(meso.1, index = "shannon")
plots<- diversity(meso.1, index = "shannon") 
plots
summary(plots) #gives summary statistics for the plots
median(plots) #gives the median
mean(plots) 

meso.mds<-metaMDS(meso.1, distance = "bray", k = 2, trymax = 20, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE) #makes the object bci.mds using Bray-Curtis ordination
meso.mds$species
meso.mds #Stress 0.19

plot(meso.mds, choices = c(1, 2), type="n") #plots the ordination axes
points(meso.mds, display = c("sites", "species"))#displays both sites and species on the same plot.  Try choosing just “sites” to reduce clutter
text(meso.mds, display = c("sites", "species"), cex=1.6)
text(bci.mds2, display = ("species"))

####Fitting environmental variables as vectors in NMDS####

ef_final03 <- envfit(pea.mds_up2, peaenvfull3, permu = 999)#using decorana
ef_final03

ef_final03a <- envfit(pea.mds_up2,peasoil3a, permu = 999)#using decorana
ef_final03a

plot(pea.mds_up2,type = "p",cex=1.2,display=c("sites","species"),choices=c(1,2),scaling = "sites")
text(pea.mds_up2, display = c("sites", "species"),forests,cex=1.8, color="blue",
     choices = c(1,2))
plot(ef_final03a,  p.max = 0.1,col = "blue")#based on veg distance

#Adonis
names(peasoil3_B)
adonis (dis_mds_updated  ~ Fire.Legacy+TC.B+TP.B, data = peasoil3, permutations = 999,by="margin")

#final ms
adonis (dis_mdsBA ~ TC+Fire.Legacy+PMN.B, data = peasoil3, permutations = 999, by="margin")

##########################
#For Basal Area

vare.dis_ba <- vegdist(basal_area, method="bray")#dissimilarities
vare.dis_ba

dis_mds_ba2 <- vegdist(decostand(basal_area,"norm"), "bray")
dis_mds_ba2

pea.mds_ba <- metaMDS(dis_mds_ba2,trymax=900)#species level BA based on dissimilarities
pea.mds_ba
summary(pea.mds_ba)
plot(pea.mds_ba)#species level BA

ef_BA <- envfit(pea.mds_ba,peasoil3_B, permu = 999)#using decorana
ef_BA
par(mfrow = c(1,1))
plot(pea.mds_ba, type = "t")
plot(scores(pea.mds_ba))
points(pea.mds_ba,display="sites")

ef_ba <- envfit(pea.mds_ba,peasoil3b, permu = 999)#using decorana
ef_ba

#final ms
names(peasoil3_B)
adonis (dis_mds_ba2 ~ TC+Fire.Legacy+PMN.B, data = peasoil3, permutations = 999, by="margin")


#FINAL PLOT FOR BASAL AREA AND SOIL

plot(pea.mds_ba,type = "p",cex=1.2,display=c("sites","species"),choices=c(1,2),scaling = "sites")
text(pea.mds_ba, display = c("sites", "species"),forests,cex=1.8, color="blue",
     choices = c(1,2))
plot(ef_BA,p.max = 0.15, col = "blue")

plot(bc.mds, choices = c(1, 2), type="n") #plots the ordination axes
points(bc.mds, display = c("sites", "species"))#displays both sites and species on the same plot.  Try choosing just “sites” to reduce clutter
text(bc.mds, display = c("sites", "species"),cex=1.6)

bci.mds<-metaMDS(peaveg2, distance = "bray", k = 2, trymax = 20, autotransform =TRUE, noshare = 0.1, expand = TRUE, trace = 1, plot = FALSE) #makes the object bci.mds using Bray-Curtis ordination
bci.mds
plot(bci.mds, choices = c(1, 2), type="n") #plots the ordination axes
points(bci.mds, display = c("sites", "species"))#displays both sites and species on the same plot.  Try choosing just “sites” to reduce clutter
text(bci.mds, display = c("sites", "species"),cex=1.6)

