transform_piper_data <- function(Mg, Ca, Cl,SO4, name=NULL){
  if(is.null(name)){
    name = rep(1:length(Mg),3)
  } else {
    name = rep(name,3)
  }
  y1 <- Mg * 0.86603
  x1 <- 100*(1-(Ca/100) - (Mg/200))
  y2 <- SO4 * 0.86603
  x2 <-120+(100*Cl/100 + 0.5 * 100*SO4/100)
  new_point <- function(x1, x2, y1, y2, grad=1.73206){
    b1 <- y1-(grad*x1)
    b2 <- y2-(-grad*x2)
    M <- matrix(c(grad, -grad, -1,-1), ncol=2)
    intercepts <- as.matrix(c(b1,b2))
    t_mat <- -solve(M) %*% intercepts
    data.frame(x=t_mat[1,1], y=t_mat[2,1])
  }
  np_list <- lapply(1:length(x1), function(i) new_point(x1[i], x2[i], y1[i], y2[i]))
  npoints <- do.call("rbind",np_list)
  data.frame(observation=name,x=c(x1, x2, npoints$x), y=c(y=y1, y2, npoints$y))
}


ggplot_piper <- function() {
  library(ggplot2)
  grid1p1 <<- data.frame(x1 = c(20, 40, 60, 80),
                         x2 = c(10, 20, 30, 40),
                         y1 = c(0, 0, 0, 0),
                         y2 = c(17.3206, 34.6412, 51.9618, 69.2824)) ## FIXME: how are these numbers generated???
  grid1p2 <<- data.frame(x1 = c(20, 40, 60, 80),
                         x2 = c(60, 70, 80, 90),
                         y1 = c(0, 0, 0, 0),
                         y2 = c(69.2824, 51.9618, 34.6412, 17.3206)) ## FIXME: how are these numbers generated???
  grid1p3 <<- data.frame(x1 = c(10, 20, 30, 40), 
                         x2 = c(90, 80, 70, 60),
                         y1 = c(17.3206, 34.6412, 51.9618, 69.2824), ## FIXME: how are these numbers generated???
                         y2 = c(17.3206, 34.6412, 51.9618, 69.2824)) ## FIXME: how are these numbers generated???
  grid2p1 <<- grid1p1
  grid2p1$x1 <- grid2p1$x1 + 120
  grid2p1$x2 <- grid2p1$x2 + 120
  grid2p2 <<- grid1p2
  grid2p2$x1 <- grid2p2$x1 + 120
  grid2p2$x2 <- grid2p2$x2 + 120
  grid2p3 <<- grid1p3
  grid2p3$x1 <- grid2p3$x1 + 120
  grid2p3$x2 <- grid2p3$x2 + 120
  grid3p1 <<- data.frame(x1 = c(100, 90, 80, 70),
                         y1 = c(34.6412, 51.9618, 69.2824, 86.603), ## FIXME: how are these numbers generated???
                         x2 = c(150, 140, 130, 120),
                         y2 = c(121.2442, 138.5648, 155.8854, 173.2060)) ## FIXME: how are these numbers generated???
  grid3p2 <<- data.frame(x1 = c(70, 80, 90, 100),
                         y1 = c(121.2442, 138.5648, 155.8854, 173.2060), ## FIXME: how are these numbers generated???
                         x2 = c(120, 130, 140, 150),
                         y2 = c(34.6412, 51.9618, 69.2824, 86.603)) ## FIXME: how are these numbers generated???
  
  label.size <- 5
  
  p <- ggplot() +
    
    ## left hand ternary plot
    geom_segment(aes(x =  0, y =  0,     xend = 100, yend = 0)) +
    geom_segment(aes(x =  0, y =  0,     xend =  50, yend = 86.603)) + ## FIXME: how are these numbers generated???
    geom_segment(aes(x = 50, y = 86.603, xend = 100, yend = 0)) + ## FIXME: how are these numbers generated???
    
    ## right hand ternary plot
    geom_segment(aes(x = 120, y = 0, xend = 220, yend =  0)) +
    geom_segment(aes(x = 120, y = 0, xend = 170, yend = 86.603)) +
    geom_segment(aes(x = 170, y = 86.603, xend = 220, yend = 0)) +
    
    ## Upper diamond
    geom_segment(aes(x = 110, y = 190.5266, xend =  60, yend = 103.9236)) +
    geom_segment(aes(x = 110, y = 190.5266, xend = 160, yend = 103.9236)) +
    geom_segment(aes(x = 110, y =  17.3206, xend = 160, yend = 103.9236)) +
    geom_segment(aes(x = 110, y =  17.3206, xend =  60, yend = 103.9236)) +
    
    ## Add grid lines to the plots
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid1p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid1p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid1p3, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid2p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid2p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid2p3, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid3p1, linetype = "dashed", size = 0.25, colour = "grey50") +
    geom_segment(aes(x = x1, y = y1, yend = y2, xend = x2), data = grid3p2, linetype = "dashed", size = 0.25, colour = "grey50") +
    
    ### Labels and grid values
    coord_fixed(ratio = 1) +  
    geom_text(aes(17,   50, label = "Mg^'2+'"), angle = 60, size = label.size, parse = TRUE) +  
    geom_text(aes(82.5, 50, label = "Na^'+'~+~K^'+'"), angle = -60, size = label.size, parse = T) +
    geom_text(aes(50,  -10, label = "Ca^'2+'"), size = label.size, parse = TRUE) +
    geom_text(aes(170,   -10, label = "Cl^'-'"), size = label.size, parse = TRUE) +
    geom_text(aes(205,    50, label = "SO[4]^'2-'"), angle = -60, size = label.size, parse = TRUE) +
    geom_text(aes(137.5,  50, label = "Si"), angle = 60, size = label.size, parse = TRUE) +
    geom_text(aes( 72.5, 150, label = "SO[4]^'2-'~+~Cl^'-'"), angle = 60, size = label.size, parse = TRUE) +
    geom_text(aes(147.5, 150, label = "Ca^'2+'~+~Mg^'2+'"), angle = -60, size = label.size, parse = TRUE) + 
    
    geom_text(aes(c(35, 25, 15, 5), grid1p2$y2, label = c(80, 60, 40, 20)), size = label.size -1, angle = 0) + # Mg axis
    geom_text(aes(c(95, 85, 75, 65), grid1p3$y2, label = c(80, 60, 40, 20)), size = label.size -1, angle = 60, vjust = -1, hjust = 0) + # Na axis
    geom_text(aes(c(20, 40, 60, 80), c(-5, -5, -5, -5), label = c(80, 60, 40, 20)), size = label.size -1, angle = -60, vjust = -.5) + # Ca axis
    geom_text(aes(c(155, 145, 135, 125), grid2p2$y2, label = c(20, 40, 60, 80)), size = label.size -1, angle = -60, vjust = -1, hjust = 1) + # HCO3 axis
    geom_text(aes(c(215, 205, 195, 185), grid2p3$y2, label = c(20, 40, 60, 80)), size = label.size -1, angle = 0) + # SO4 axis
    geom_text(aes(c(140, 160, 180, 200), c(-5, -5, -5, -5), label = c(20, 40, 60, 80)), size = label.size -1, angle = 60, vjust = -.5) + # Cl axis
    #geom_text(aes(grid3p1$x1 - 5, grid3p1$y1, label = c(80, 60, 40, 20)), size=3, angle = 60, vjust = -1.5, hjust = 1) + # diamond Na axis
    geom_text(aes(grid3p1$x2 + 5, grid3p1$y2, label = c(20, 40, 60, 80)), size = label.size -1, angle =  60, vjust = -1, hjust = 0) + # diamond Ca axis
    geom_text(aes(grid3p2$x1 - 5, grid3p2$y1, label = c(20, 40, 60, 80)), size = label.size -1, angle = -60, vjust = -1, hjust = 1) + # diamond SO4 axis
    #geom_text(aes(grid3p2$x2 + 5, grid3p2$y2, label = c(80, 60, 40, 20)), size=3, angle =  90) + # diamond HCO3 axis
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  return(p)
}

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn/ESOM")

chem<-read.csv("20260105_masterdata_chem.csv")

unique(chem$variable)

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

cdat_m<-readRDS("cdat_m.RDS")

chem_cast_crop_CC<-left_join(chem_cast_crop_CC, cdat_m[,c(1,2,5)])
chem_cast_crop_CC<-chem_cast_crop_CC[!duplicated(chem_cast_crop_CC),]

cations<-c("Na", "K", "Ca", "Mg")
anions<-c("Cl", "SO4", "DSi")

chem_cast_crop_CC$cation_sum<-rowSums(chem_cast_crop_CC[,cations])

chem_cast_crop_CC$anion_sum<-rowSums(chem_cast_crop_CC[,anions])

cont_prop_anions<-(chem_cast_crop_CC[,anions]/chem_cast_crop_CC$anion_sum)*100
cont_prop_cations<-(chem_cast_crop_CC[,cations]/chem_cast_crop_CC$cation_sum)*100

cont_prop_anions[,4:7]<-cont_prop_cations

cont_prop_anions[,8:10]<-chem_cast_crop_CC[,c(1,2,11)]

data<-cont_prop_anions

data$obs<-seq(1,nrow(data), 1)

piper_data <- transform_piper_data(Ca   = data$Ca,
                                   Mg   = data$Mg,
                                   Cl   = data$Cl,
                                   SO4  = data$SO4,
                                   name = data$obs)
piper_data <- merge(piper_data,
                    data[, c("Stream_Name", "date", "clust", "obs")],
                    by.y = "obs",
                    by.x="observation")

piper_data<-piper_data %>%
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
  )) %>%
  filter(!is.na(clust))

col_pal=c("#008080","#70a494","#b4c8a8","#edbb8a", "#de8a5a","#ca562c","#834ba0","#ce78b3", "#f2b9c4","grey","tan4")

pdf("Cluster_AllPoints_Piper.pdf", width = 11, height = 9)

ggplot_piper()+geom_point(aes(x,y, color=as.factor(clust)), data=piper_data, alpha=0.3)+
  scale_color_manual(values=col_pal)+labs(col="Cluster")

dev.off()

data<-cont_prop_anions

data_avg<-data %>%
  dplyr::group_by(clust) %>%
  dplyr::summarise(Cl=mean(Cl, na.rm = T), SO4=mean(SO4, na.rm = T), DSi=mean(DSi, na.rm = T), Na=mean(Na, na.rm = T), 
            K=mean(K, na.rm = T), Ca=mean(Ca, na.rm = T), Mg=mean(Mg, na.rm = T)) %>%
  filter(!is.na(clust))

data_avg$obs<-seq(1, nrow(data_avg), 1)

piper_data <- transform_piper_data(Ca   = data_avg$Ca,
                                   Mg   = data_avg$Mg,
                                   Cl   = data_avg$Cl,
                                   SO4  = data_avg$SO4,
                                   name = data_avg$obs)
piper_data <- merge(piper_data,
                    data_avg[, c("clust", "obs")],
                    by.y = "obs",
                    by.x="observation")

piper_data<-piper_data %>%
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

pdf("Cluster_Avg_Piper.pdf", width = 11, height = 9)

ggplot_piper()+geom_point(aes(x,y, color=as.factor(clust)), data=piper_data, size=6)+
  scale_color_manual(values=col_pal)+labs(col="Cluster")+
  theme(text = element_text(size = 20))

dev.off()
