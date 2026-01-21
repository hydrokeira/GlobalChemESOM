####Cluster ESOM####
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

res<-readRDS("SiESOM_200epochs_01182026_75_120.RDS")

#Make and plot the UMatrix
umatrix <-  umatrixForEsom(res$Weights, Lines=res$Lines, 
                           Columns=res$Columns, Toroid=res$Toroid)

#add unique ID for BMU
BestMatches <- as.data.frame(res$BestMatches)
#row of BM -1 times number of columns (82) plus column of BM (not sure why its like this..)
ids <- ((res$BestMatches[,2]-1)*length(unique(res$BestMatches[,3])))+res$BestMatches[,3]
#add unique ID to BMU
BestMatches<-cbind(BestMatches,ids)
names(BestMatches) <- c("obs", "row", "col", "id")

weights_clust<-res$Weights

#kmeans using 6 clusters
set.seed(123)

p1<-fviz_nbclust(weights_clust, kmeans, method="wss", k.max = 20)

p1

p2<-fviz_nbclust(weights_clust, kmeans, method="silhouette", k.max = 20)

p2

pdf("Cluster_Metrics_Kmeans_10182025.pdf", width = 7, height = 8)

ggarrange(p2, p1, nrow=2, align = "v")

dev.off()

set.seed(123)
kmeans_cluster<-kmeans(weights_clust, iter.max=50, nstart=50, centers = 11)

clusts<-kmeans_cluster$cluster

kmeans_clusts_mat<-matrix(clusts, nrow = 75, ncol=120, byrow=TRUE)

kmeans_clusts_mat_melt<-melt(kmeans_clusts_mat)

colnames(kmeans_clusts_mat_melt)<-c("row", "col", "clust")

BestMatches_cluster<-merge(BestMatches, kmeans_clusts_mat_melt, by=c("row", "col"))

#add clusters to data and plot composition of each cluster
cont_prop$obs <- seq(1,nrow(cont_prop),1)
cdat <- inner_join(cont_prop,BestMatches_cluster[,c("obs", "clust")], by=c("obs"))

dis <- dist(cdat[,solutes_2])^2 #squared euclidean distance
#dis <- vegdist(cdat[,param], method="euclidean")
sil <- silhouette(cdat$clust, dis)
# 
# sil$sil_width

#windows()
#plot(sil)

pdf("Sil_Width_Jan2024.pdf", width = 8, height = 4, family = "Times")

fviz_silhouette(sil)+
  labs(title="", fill="Cluster", col="Cluster", y="Silhouette width")
scale_color_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                            "4" = colors[4], "5"=colors[5], "6"=colors[6]))+
  scale_fill_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                             "4" = colors[4], "5"=colors[5], "6"=colors[6]))

dev.off()

BestMatches_cluster <- BestMatches_cluster %>%
  mutate(clust=case_when(
    clust==5~1,
    clust==10~2,
    clust==8~3,
    clust==11~4,
    clust==1~5,
    clust==7~6,
    clust==4~7,
    clust==2~8,
    clust==3~9,
    clust==9~10,
    clust==6~11
  ))

col_pal=c("#008080","#70a494","#b4c8a8","#edbb8a", "#de8a5a","#ca562c","#834ba0","#ce78b3", "#f2b9c4","grey","tan4")

#f9ddda,#f2b9c4,#e597b9,#ce78b3,#ad5fad,#834ba0,#573b88

pdf("ESOM_Matrix_75_120_11Clusters_01212026_smallpoints.pdf", width = 8, height = 6)

#plot umatrix with clusters on top
plotMatrix(umatrix, as.matrix(BestMatches_cluster[,c("obs","row","col")]), 
           Toroid=F, BmSize=3, DrawLegend=T, Clean=T, 
           Cls=as.numeric(BestMatches_cluster$clust), ClsColors = col_pal)


dev.off()

id.vars<-c("Stream_Name", "date", "check", "obs", "clust")

cdat_m <- reshape2::melt(cdat, id.vars=id.vars)
labs <- c(paste("Cluster", seq(1,11,1)))

cdat_m$variable <- as.character(cdat_m$variable)
cdat_m$variable<-factor(cdat_m$variable, levels = c("Ca", "Mg", "Na", "K", "N", "SO4", "Cl", "DSi"))
cdat_m$value<-as.numeric(cdat_m$value)

cdat_m <- cdat_m %>%
  mutate(clust=case_when(
    clust==5~1,
    clust==10~2,
    clust==8~3,
    clust==11~4,
    clust==1~5,
    clust==7~6,
    clust==4~7,
    clust==2~8,
    clust==3~9,
    clust==9~10,
    clust==6~11
  ))

pdf("ESOM_Cluster_ClusterConcentrationBoxplots_10162025.pdf", width = 10, height = 12)

ggplot(data=cdat_m, aes(x=variable, y=value, fill=as.factor(clust)))+
  geom_violin(position=position_dodge(width=1.2), scale = "width")+
  labs(x="Solute", y="Proportion of concentration")+
  theme_bw() + 
  #scale_fill_manual(values = carto_pal(n=11, "Vivid"))+
  scale_fill_manual(values = col_pal)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust = 1),legend.position="none", text = element_text(size=20))+
  facet_wrap(~clust, nrow=4)

dev.off()

clust_stream <- cdat %>%
  group_by(clust) %>%
  distinct(Stream_Name) %>%
  tally()

stream_clust <- cdat %>%
  group_by(Stream_Name) %>%
  distinct(clust)

stream_clust_count<-
  stream_clust %>%
  group_by(Stream_Name) %>%
  tally()

clust_stream <- clust_stream %>%
  mutate(clust=case_when(
    clust==5~1,
    clust==10~2,
    clust==8~3,
    clust==11~4,
    clust==1~5,
    clust==7~6,
    clust==4~7,
    clust==2~8,
    clust==3~9,
    clust==9~10,
    clust==6~11
  ))

p2<-ggplot(stream_clust_count, aes(x=n))+geom_bar(fill="black")+
  labs(x="Number of Cluster Memberships", y="Number of Sites")+
  scale_x_continuous(labels=seq(1,9,1), breaks = seq(1,9,1))+
  theme_classic()+theme(text = element_text(size = 20))

p1<-ggplot(clust_stream, aes(x=clust, y=n))+geom_bar(fill="black", stat = "identity")+
  labs(x="Cluster", y="Number of Sites")+theme_classic()+
  scale_x_continuous(labels=seq(1,13,1), breaks = seq(1,13,1))+
  theme(text = element_text(size = 20))

pdf("ClusterStats_01212026.pdf", width = 6, height = 8)

ggarrange(p1, p2, nrow = 2)

dev.off()
