#Applied Soil Ecology Manuscript Figures

#Soil Mesofauna Diversity and Abundance

library(ggplot2)
library(patchwork)
library(ggcorrplot)
library(ggpubr)

meso_cap2<-read.csv(file.choose())#20210502_FAUNA_NMDS_ASE
str(meso_cap2)
names(meso_cap2)
meso_cap2.1 = as.data.frame(sapply(meso_cap2[3:27], as.numeric))
str(meso_cap2.1)

H<-diversity(meso_cap2.1, index = "shannon")
H
S <- specnumber(meso_cap2.1) ## rowSums(BCI > 0) does the same...
S
J <- H/log(S)#Pielou's evenness J = H/log(S)
J
Plot<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
Plot
comm_data<-cbind(data.frame(meso_cap2$Season,meso_cap2$Plot,Plot,meso_cap2$Total_Abundance,H,S,J))
names(comm_data)
comm_data

#exporting comm-data as csv /Users/bbomfim/Dropbox/R/Inkotte_Chap3
write.csv(comm_data,"C:Users\\bbomfim\\Dropbox\\R\\Inkotte_Chap3\\Community_Data_Cap2.csv", row.names = FALSE)

p_div<-ggplot(comm_data, aes(x=reorder(meso_cap2.Season, H),y=H, fill = meso_cap2.Season)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +  #this flips the co-ordinates so your x axis becomes your y and vice versa
  labs(y="Shannon's diversity", x="") +#scale_color_discrete(palette="Dark2")
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(panel.grid.minor=element_blank(), # get rid of background grid
        panel.grid.major=element_blank(),
        axis.title=element_text(size=20),
        axis.text=element_text(size=19),legend.position="none")+scale_fill_grey(start = 0.2,
                                                                                end = 0.8)# use theme() to change font sizes
p_div2<-p_div + stat_compare_means(method = "t.test",aes(label = ..p.signif..), 
                                                                                label.x = 1.98, label.y = 1.78,cex=10)
p_div2
p_div3<-p_div + stat_compare_means( 
                                  label.x = 2, label.y = 1.9,cex=6)+ stat_compare_means(method = "t.test",aes(label = ..p.signif..), 
                                                                                        label.x = 1.98, label.y = 1.78,cex=10)
p_div3
p_abd<-ggplot(comm_data, aes(x=reorder(meso_cap2.Season, -meso_cap2.Total_Abundance),y=meso_cap2.Total_Abundance, fill = meso_cap2.Season)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +  #this flips the co-ordinates so your x axis becomes your y and vice versa
  labs(y="Total abundance", x="Season") +#scale_color_discrete(palette="Dark2")
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(panel.grid.minor=element_blank(), # get rid of background grid
        panel.grid.major=element_blank(),
        axis.title=element_text(size=20),
        axis.text=element_text(size=19),legend.position="none")+scale_fill_grey(start = 0.2,
                                                                                end = 0.8)# use theme() to change font sizes
p_abd
p_abd2<-p_abd +  stat_compare_means(method = "t.test",aes(label = ..p.signif..),label.x = 1, label.y = 600,cex=8)
p_abd2
p_abd3<-p_abd + stat_compare_means( 
                                   label.x = 2, label.y = 600,cex=6)+ stat_compare_means(method = "t.test",aes(label = ..p.signif..), 
                                                                                         label.x = 2, label.y = 480,cex=8)
p_abd3

p_even<-ggplot(comm_data, aes(x=reorder(meso_cap2.Season, J),y=J, fill = meso_cap2.Season)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +  #this flips the co-ordinates so your x axis becomes your y and vice versa
  labs(y="Pielou's evenness", x="Season") +#scale_color_discrete(palette="Dark2")
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(panel.grid.minor=element_blank(), # get rid of background grid
        panel.grid.major=element_blank(),
        axis.title=element_text(size=20),
        axis.text=element_text(size=19),legend.position="none")+scale_fill_grey(start = 0.2,
                                                                                end = 0.8)# use theme() to change font sizes
p_even
p_even2<-p_even + stat_compare_means(method = "t.test",aes(label = ..p.signif..), label.x = 2.02, label.y = 0.83,cex=10)
p_even2
p_even3<-p_even + stat_compare_means( 
                                     label.x = 2.2, label.y = 0.77,cex=6)+ stat_compare_means(method = "t.test",aes(label = ..p.signif..), 
                                                                                              label.x = 2, label.y = 0.85,cex=10)
p_even3

p_sp<-ggplot(comm_data, aes(x=reorder(meso_cap2.Season, -S),y=S, fill = meso_cap2.Season)) +
  geom_boxplot(alpha=0.8) +
  coord_flip() +  #this flips the co-ordinates so your x axis becomes your y and vice versa
  labs(y="Number of groups", x="") +#scale_color_discrete(palette="Dark2")
  theme_bw() + # ggplot2 has a few theme options, I like minimal and classic
  theme(panel.grid.minor=element_blank(), # get rid of background grid
        panel.grid.major=element_blank(),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size=19),axis.text.x=element_text(hjust=0.9,size=19),legend.position="none")+scale_fill_grey(start = 0.2,
                                                                                end = 0.8)# use theme() to change font sizesp_sp
p_sp
p_sp2<-p_sp + stat_compare_means(method = "t.test",aes(label = ..p.signif..), label.x = 1.12, label.y = 16,cex=8)
p_sp2
p_sp3<-p_sp + stat_compare_means( 
                                 label.x = 1.22, label.y = 16,cex=6)+ stat_compare_means(method = "t.test",aes(label = ..p.signif..), 
                                                                                        label.x = 1.1, label.y = 14,cex=8)
p_sp3

#Final plot
plot_meso_final<-(p_abd2+p_div2)/(p_even2+p_sp2)+
plot_meso_final

ggsave(filename = "Fig_Meso_Cap2_v5.png",
       plot = plot_meso_final, width = 18, height = 13, units = 'cm',
       scale = 2, dpi = 1000)


#### Correlation plot####

soilcor<-read.csv(file.choose())#20210503_Soil
soilcor$Total_abundance=as.numeric(soilcor$Total_abundance)
attach(soilcor)
head(soilcor)
names(soilcor)
names(comm_data)
soilcor$S<-comm_data$S
soilcor$J<-comm_data$J

Figdry <- soilcor
Figdry2<-Figdry[,c(3:12)]
names(Figdry2)
names(Figdry2)[names(Figdry2) == "Total_abundance"] <- "Total abundance"
names(Figdry2)[names(Figdry2) == "Shannon_diversity_h"] <- "Shannon's diversity"
names(Figdry2)[names(Figdry2) == "Total_litterfall_g_m2_month"] <- "Montly total litterfall"
names(Figdry2)[names(Figdry2) == "Leaf_fall_g_m2_month"] <- "Monthly leaf fall"
names(Figdry2)[names(Figdry2) == "Wood_fall_g_m2_month"] <- "Monthly wood fall"
names(Figdry2)[names(Figdry2) == "Mean_monthly_rainfall_mm"] <- "Mean monthly rainfall"
names(Figdry2)[names(Figdry2) == "Leaf_k_at_720"] <- "Leaf decomposition rate"
names(Figdry2)[names(Figdry2) == "Wood_k_at_720"] <- "Wood decomposition rate"
names(Figdry2)[names(Figdry2) == "S"] <- "Number of groups"
names(Figdry2)[names(Figdry2) == "J"] <- "Pielou's evenness"

corrdry <- round(cor(Figdry2,method="pearson"), 2)

p.matdry <- cor_pmat(Figdry2)

ggcorrplot(corrdry)
ggcorrplot(corrdry, hc.order = TRUE, type = "lower",hc.method = "ward.D2",
           outline.col = "white", p.mat = p.matdry,method="circle",ggtheme=ggplot2::theme_bw())

Fig3_corr<-ggcorrplot(corrdry, hc.order = TRUE, type = "lower",hc.method = "ward.D2",sig.level = 0.05,
                    outline.col = "white", p.mat = p.matdry,method="square",ggtheme=ggplot2::theme_minimal(),show.legend=TRUE, 
                    legend.title="Pearson's r", lab=TRUE, lab_size=6, tl.cex=26,
                    colors = c("#ABA0A0", "white", "#ffa600",pch.cex=20,nbreaks = 8,legend.text.cex=26))+font("legend.text",size=16)+font("legend.title", size=18)#+theme(axis.text.x = element_text(margin=margin(-2,0,0,0)),axis.text.y = element_text(margin=margin(0,-2,0,0)))
Fig3_corr

#Saving figure in high res
ggsave(filename = "Fig3_Final.png",
       plot = Fig3_corr, width = 16, height = 18, units = 'cm',
       scale = 2, dpi = 1000)

##END###
