##PCA soil, litterfall, rainfall##

library(devtools)
install_github("vqv/ggbiplot")
library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(ggbiplot)


#uploading data
PCA_ms<-read.csv(file.choose())#20210503_Cap2_Env_Matrix
attach(PCA_ms)
names(PCA_ms)
str(PCA_ms)#n=12 of 23 variables
summary(PCA_ms)
PCA_ms$K=as.numeric(PCA_ms$K)

#changing column names
names(PCA_ms)[names(PCA_ms) == "t"] <- "CECe"
names(PCA_ms)[names(PCA_ms) == "T"] <- "CECt"
names(PCA_ms)[names(PCA_ms) == "Mean_Richness"] <- "Richness"
names(PCA_ms)[names(PCA_ms) == "Annual_Abundance"] <- "Abundance"
names(PCA_ms)[names(PCA_ms) == "Mean_Diversity"] <- "Diversity"
names(PCA_ms)[names(PCA_ms) == "Mean_Evenness"] <- "Evenness"
names(PCA_ms)[names(PCA_ms) == "P.Rem"] <- "Prem"
names(PCA_ms)[names(PCA_ms) == "P"] <- "Pav"
names(PCA_ms)[names(PCA_ms) == "C_stock"] <- "SOCstock"
names(PCA_ms)[names(PCA_ms) == "H.Al"] <- "H+Al"
names(PCA_ms)[names(PCA_ms) == "OM"] <- "SOM"

#checking new names
names(PCA_ms)

#data frame
finalpca<-PCA_ms[,c(2:5,10:12,16,19:21,23)]
names(finalpca)

#running prcomp function
mtcars.pca <- prcomp(finalpca, center = TRUE,scale. = TRUE)#standardized data
summary(mtcars.pca)
str(mtcars.pca)
mtcars.pca$rotation

#Inspecting results
mtcars.pca$sdev
#The values of each sample in terms of the principal components ($x)
mtcars.pca$x
#The relationship (correlation or anticorrelation, etc) between the initial variables and the principal components ($rotation)
mtcars.pca$rotation

head(finalpca)
plots<-PCA_ms[,1]
plots

#Plotting biplot
pca_final<-ggbiplot(mtcars.pca, varname.size=5.5,varname.adjust=1.5,color="darkgray",labels.size=4)+geom_point(shape=21,size=2,fill="white",stroke=0.8,color="black")+
  theme_classic()+ theme(plot.title = element_text(size=18,hjust = 0.5))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+xlim(-2,2.2)+#geom_point(aes(x=plots,y=plots),shape = 21,size = 5.5)+
  theme(axis.title = element_text (size = 21,face="bold"), 
        axis.text = element_text (size = 20,face="bold"))+theme(plot.title = element_text(size=18,hjust = 0.5))+
  labs(x="PC1 (47.4% explained variance)",y="PC2 (19% explained variance)",title="PCA soil and litter")#+ggrepel::geom_label_repel(aes(label=plots),box.padding = 0.8,cex = 5.5, direction = "both", segment.size = 0.25)

pca_final

##END##