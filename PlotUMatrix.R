#clean data run ESOM
#plot ESOM data

require(reshape2)
require(tidyverse)
library(kohonen)
library(Umatrix)
library(gridExtra)
library(factoextra)
library(ggpubr)
library(cluster)
library(rcartocolor)
library(reshape2)
library(dplyr)
library(cetcolor)
library(lemon)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn/ESOM")

chem<-read.csv("20260105_masterdata_chem.csv")

solutes<-c("Ca", "Cl", "DSi", "K", "Mg", "Na", "NO3", "NOx", "SO4")

solutes_2<-c("Ca", "Cl", "DSi", "K", "Mg", "Na", "N", "SO4")

chem<-subset(chem, chem$variable %in% solutes)

chem_cast<-dcast(chem, Stream_Name+date~variable, value.var = "value", fun.aggregate = mean)

chem_cast[,3:11]<-lapply(chem_cast[,3:11], as.numeric)

na_count <- data.frame(sapply(chem_cast, function(y) sum(length(which(is.na(y))))) )

chem_cast$N<-rowMeans(chem_cast[,c(9,10)], na.rm = T)

chem_cast<-chem_cast[,c(1:8,11,12)]

chem_cast_crop_CC<-chem_cast[complete.cases(chem_cast),]

all_sites<-chem_cast_crop_CC %>%
  group_by(Stream_Name) %>%
  summarise(n_obs=n())

good_sites<-chem_cast_crop_CC %>%
  group_by(Stream_Name) %>%
  summarise(n_obs=n()) %>%
  filter(n_obs > 9) %>%
  filter(!str_detect(Stream_Name, "FDUP"))

chem_cast_crop_CC <- chem_cast_crop_CC %>%
  filter(Stream_Name %in% good_sites$Stream_Name)

chem_cast_crop_CC$sum<-rowSums(chem_cast_crop_CC[,3:10])

cont_prop<-chem_cast_crop_CC[,3:10]/chem_cast_crop_CC$sum

cont_prop<-cbind(chem_cast_crop_CC[c(1,2)], cont_prop)

cont_prop$check<-rowSums(cont_prop[3:10])

cont_prop<-cont_prop[complete.cases(cont_prop),]

#Normalizing data - needs to save mean and variance for later
z.means <- apply(cont_prop[,solutes_2],2,mean, na.rm=TRUE)
z.vars <- apply(cont_prop[,solutes_2],2, sd, na.rm=TRUE)
#write.csv (rbind(z.means,z.vars),"Data means and vars.csv",row.names=FALSE)

X1 <- scale(cont_prop[,solutes_2], center=TRUE, scale = TRUE) #These are the normalized input data for the ESOM

#only run this line on version that is run on the server
#res <-  esomTrain(X1, InitMethod = "norm_mean_std2", Epochs=100, Key = 1:nrow(X1), Lines = 75, Columns = 120)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn/ESOM")

#only run this line on version that is run on the server
#saveRDS(res, "SiESOM_11242025_75_120.RDS")

res<-readRDS("SiESOM_200epochs_01182026_75_120.RDS")

#Make and plot the UMatrix
umatrix <-  umatrixForEsom(res$Weights, Lines=res$Lines, 
                           Columns=res$Columns, Toroid=res$Toroid)

pdf("ESOM_Matrix_Final_75_120_01182026.pdf", width = 11, height = 8.5)

plotMatrix(umatrix, Toroid=FALSE, BmSize = 2,
           DrawLegend = T, Clean = TRUE)+theme(text = element_text(size=20))

dev.off()

pdf("ESOM_Matrix_Torroid_Final_75_120_01182026.pdf", width = 11, height = 8.5)

plotMatrix(umatrix, Toroid=TRUE, BmSize = 2,
           DrawLegend = T, Clean = TRUE)+theme(text = element_text(size=20))

dev.off()


#add unique ID for BMU
BestMatches <- as.data.frame(res$BestMatches)
#row of BM -1 times number of columns (82) plus column of BM (not sure why its like this..)
ids <- ((res$BestMatches[,2]-1)*length(unique(res$BestMatches[,3])))+res$BestMatches[,3]
#add unique ID to BMU
BestMatches<-cbind(BestMatches,ids)
names(BestMatches) <- c("obs", "row", "col", "id")

#Convert data to unnormalized results
#str(res)
Wts.raw <- res$Weights
Wts.unnorm <- (z.vars)*t(Wts.raw)+z.means
Wts.unnorm <- t(Wts.unnorm)

colnames(Wts.unnorm) <- colnames(cont_prop[,solutes_2]) #move the column names over

Wts.unnorm<-Wts.unnorm[,c("Ca", "Mg", "Na", "K", "N", "SO4", "Cl", "DSi")]

#plot unnormalized compositional weights for each solute
pplots <- list()
x<-1
for(i in 1:dim(Wts.unnorm)[2]){
  cname <- colnames(Wts.unnorm)[i]
  test <- matrix(Wts.unnorm[,i], nrow = 75, ncol=120, byrow=TRUE)
  pplots[[x]] <- plotMatrix(test, Toroid=FALSE, Clean = TRUE, Title=cname, DrawLegend=F, ColorStyle = "Pmatrix")+
    theme(text = element_text(size = 20))
  x <- x+1
}

plot_legend<-plotMatrix(test, Toroid=FALSE, Clean = TRUE, Title=cname, DrawLegend=T, ColorStyle = "Pmatrix")+
  theme(legend.position = "bottom", text = element_text(size = 20), legend.text = element_text(angle = 45, hjust = 1))+labs(fill="Compositional Solute Concentration")+
  guides(fill=guide_colorbar(title.position = "top"))

plot_legend

leg<-g_legend(plot_legend)

hlay <- rbind(c(1,1,1,2,2,2),
              c(1,1,1,2,2,2),
              c(4,4,4,5,5,5),
              c(4,4,4,5,5,5),
              c(6,6,6,7,7,7),
              c(6,6,6,7,7,7),
              c(8,8,8,9,9,9),
              c(8,8,8,9,9,9),
              c(3,3,3,3,3,3))

grid.arrange(pplots[[1]],pplots[[2]],leg, pplots[[3]],pplots[[4]],pplots[[5]],
             pplots[[6]],pplots[[7]],pplots[[8]], layout_matrix=hlay)

pdf("ESOM_Final_CompSolutePlots_01182026.pdf", width = 8, height = 10)

grid.arrange(pplots[[1]],pplots[[2]],leg, pplots[[3]],pplots[[4]],pplots[[5]],
             pplots[[6]],pplots[[7]],pplots[[8]], layout_matrix=hlay)

dev.off()

cont_prop$obs <- seq(1,nrow(cont_prop),1)
cdat2<-inner_join(cont_prop, BestMatches, by="obs")
cdat2$date<-as.Date(cdat2$date)
cdat2<-cdat2[complete.cases(cdat2$date),]

cdat2<-cdat2 %>%
  group_by(id) %>%
  mutate(mean_year=mean(year(date)), mean_month=mean(month(date)))

cdat2 <- cdat2 %>%
  mutate(mean_year_group = round(mean_year/10)*10)

cdat2$year_cluster<-group_indices(cdat2, mean_year_group)

colors_months<-cet_pal(12, "c4s")

colors_years<-scales::viridis_pal()(7)

p1<-plotMatrix(umatrix, as.matrix(cdat2[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=T, Clean=T, 
               Cls=round(cdat2$mean_month), ClsColors = colors_months)+ggtitle("Nodes Colored by Mean Month")+
  theme(text = element_text(size = 20))

p2<-plotMatrix(umatrix, as.matrix(cdat2[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=T, Clean=T, 
               Cls=round(cdat2$year_cluster), ClsColors = colors_years)+ggtitle("Nodes Colored by Mean Decade")+
  labs(color="Mean Decade")+theme(text = element_text(size = 20))+
  scale_color_manual(breaks = c(1,2,3,4,5,6,7), values = colors_years,
                     labels=c("1960", "1970","1980","1990", "2000", "2010", "2020"))

pdf("Temporal_Variation_Nodes_FinalESOM_01182026.pdf", width = 10, height = 12)

ggarrange(p1, p2, nrow = 2, align = "v")

dev.off()
