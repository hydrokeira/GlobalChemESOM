setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn/ESOM")

chem<-read.csv("20260105_masterdata_chem.csv")

solutes<-c("Ca", "Cl", "DSi", "K", "Mg", "Na", "NO3", "NOx", "SO4")

solutes_2<-c("Ca", "Cl", "DSi", "K", "Mg", "Na", "N", "SO4")

chem<-subset(chem, chem$variable %in% solutes)

chem_cast<-dcast(chem, Stream_Name+date~variable, value.var = "value", fun.aggregate = mean)

chem_cast[,3:11]<-lapply(chem_cast[,3:11], as.numeric)

chem_cast$N<-rowMeans(chem_cast[,c(9,10)], na.rm = T)

chem_cast<-chem_cast[,c(1:8,11,12)]

chem_cast_crop_CC<-chem_cast[complete.cases(chem_cast),]

all_sites<-chem_cast_crop_CC %>%
  group_by(Stream_Name) %>%
  summarise(n_obs=n())

good_sites<-chem_cast_crop_CC %>%
  group_by(Stream_Name) %>%
  summarise(n_obs=n()) %>%
  filter(n_obs > 9)

sites_100<-subset(good_sites, good_sites$n_obs > 99)

cdat<-readRDS("cdat.RDS")

stream_clust <- cdat %>%
  group_by(Stream_Name) %>%
  distinct(clust)

clust_stream <- cdat %>%
  group_by(clust) %>%
  distinct(Stream_Name) %>%
  tally()

stream_clust_count<-
  stream_clust %>%
  group_by(Stream_Name) %>%
  tally()

stream_clust_shifting<-subset(stream_clust_count, stream_clust_count$Stream_Name %in% sites_100$Stream_Name)

shifting_sites<-unique(stream_clust_shifting$Stream_Name)

cdat_shifting<-subset(cdat, cdat$Stream_Name %in% shifting_sites)

clust_count<-cdat_shifting %>%
  group_by(Stream_Name, clust) %>%
  count()

total_count<-cdat_shifting %>%
  group_by(Stream_Name) %>%
  count()

clust_count<-merge(clust_count, total_count, by="Stream_Name")

clust_count$prop<-(clust_count$n.x/clust_count$n.y)*100

streams_mode<-clust_count %>%
  group_by(Stream_Name) %>%
  slice_max(prop) %>%
  filter(!duplicated(Stream_Name))

streams_mode<-streams_mode[,c(1,2)]
colnames(streams_mode)[2]<-"Modal_Cluster"

cdat_shifting$date<-as.Date(cdat_shifting$date)

cdat_shifting<-merge(cdat_shifting, streams_mode, by="Stream_Name")

cdat_shifting<-left_join(cdat_shifting, stream_clust_count)

streams_mode2<-subset(streams_mode, streams_mode$Stream_Name %in% cdat_shifting$Stream_Name)

cdat_shifting <- cdat_shifting %>%
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

cdat_shifting <- cdat_shifting %>%
  mutate(Modal_Cluster=case_when(
    Modal_Cluster==5~1,
    Modal_Cluster==10~2,
    Modal_Cluster==8~3,
    Modal_Cluster==11~4,
    Modal_Cluster==1~5,
    Modal_Cluster==7~6,
    Modal_Cluster==4~7,
    Modal_Cluster==2~8,
    Modal_Cluster==3~9,
    Modal_Cluster==9~10,
    Modal_Cluster==6~11
  ))

col_pal=c("#008080","#70a494","#b4c8a8","#edbb8a", "#de8a5a","#ca562c","#834ba0","#ce78b3", "#f2b9c4","grey","tan4")

pdf("ESOM_Cluster_Shifting.pdf", width = 8, height = 12)

p1<-cdat_shifting %>%
  ggplot(aes(date, Stream_Name))+geom_point(aes(col=as.factor(clust)), alpha=0.5, size=0.8)+
  scale_color_manual(values = col_pal)+
  scale_y_discrete(limits=unique(rev(cdat_shifting$Stream_Name[order(cdat_shifting$Modal_Cluster)])))+
  theme_classic()+
  theme(axis.text.y = element_blank(),  text = element_text(size=20))+labs(y="", col="Cluster", x="Date")

dev.off()

to_binary_change <- function(x) {
  c(0, as.integer(diff(x) != 0))
}

cdat_shifting<-cdat_shifting %>%
  dplyr::group_by(Stream_Name) %>%
  mutate(binary_cluster=to_binary_change(clust))

p2<-cdat_shifting %>%
  ggplot(aes(date, Stream_Name))+geom_point(aes(col=as.factor(binary_cluster)), alpha=0.5, size=0.8)+
  scale_color_manual(values = c("grey70", "black"))+
  scale_y_discrete(limits=unique(rev(cdat_shifting$Stream_Name[order(cdat_shifting$Modal_Cluster)])))+
  theme_classic()+
  theme(axis.text.y = element_blank(), text = element_text(size=20))+labs(y="", col="Cluster", x="Date")

pdf("AllShifts_Binary.pdf", width = 14, height = 12)

ggarrange(p1, p2)

dev.off()

cdat_shifting_metrics <- cdat_shifting %>%
  dplyr::group_by(Stream_Name) %>%
  arrange(date, .by_group = TRUE) %>%
  summarise(
    n_changes=sum(clust != dplyr::lag(clust), na.rm=T),
    
    prop_modal=mean(clust == Modal_Cluster, na.rm = TRUE),
    
    n_clusters=n_distinct(clust),
    # binary series (temporary handle)
    b = list(binary_cluster),
    
    # inter-event times
    iet = list(diff(which(binary_cluster == 1))),
    
    # regularity of change (cyclicity)
    cv_iet = ifelse(
      length(iet[[1]]) > 1,
      sd(iet[[1]]) / mean(iet[[1]]),
      NA_real_
    ),
    
    # autocorrelation strength
    max_acf = max(abs(acf(binary_cluster, plot = FALSE)$acf[-1])),
    
    trend=cor(year(date), as.numeric(factor(clust)), method = "spearman"),
    
    # progressive vs stabilizing trend in change frequency
    trend_change = cor(
      seq_along(binary_cluster),
      binary_cluster,
      use = "complete.obs"
    ),
    
    .groups = "drop"
  )

col_pal=c("#008080","#70a494","#b4c8a8","#edbb8a", "#de8a5a","#ca562c","#834ba0","#ce78b3", "#f2b9c4","grey","tan4")

cdat_shifting_metrics<-cdat_shifting_metrics %>%
  mutate(shifting_classification=case_when(
    n_changes < 10 & prop_modal > 0.90 ~ "stable/mostly stable",
    cv_iet < 1.5 & abs(trend) < 0.2 & n_clusters <= 3 ~ "cyclical",
    cv_iet < 1.5 & abs(trend) < 0.2 & n_clusters > 3 ~ "non-specific cyclical",
    abs(trend) >= 0.2 & n_changes > 10 ~ "progressive",
    #abs(trend_change) > 0.15 & n_changes > 10 & abs(trend) < 0.2 ~ "stabilizing/destabilizing",
    .default="non-patterned change"
    ))

table(cdat_shifting_metrics$shifting_classification)

pdf("Shifts_Class_Facet.pdf", width = 18, height = 12)

cdat_shifting %>%
  left_join(cdat_shifting_metrics) %>%
  mutate(shifting_classification = factor(shifting_classification, 
                                          levels=c("stable/mostly stable", "progressive",
                                                   "cyclical","non-specific cyclical",
                                                   "non-patterned change"))) %>%
  #filter(shifting_classification=="progressive") %>%
  dplyr::group_by(Stream_Name) %>%
  mutate(
    Modal_Cluster_stream = unique(Modal_Cluster)  # or unique()
  ) %>%
  ungroup() %>%
  mutate(
    Stream_Name = fct_reorder(Stream_Name, Modal_Cluster_stream, .desc = TRUE)
  ) %>%
  ggplot(aes(date, Stream_Name))+geom_point(aes(col=as.factor(clust)), alpha=0.8, size=2)+
  scale_color_manual(values = col_pal)+
  theme_classic()+
  theme(axis.text.y = element_blank(), text = element_text(size = 20))+labs(y="", col="Cluster", x="Date")+
  facet_wrap(~shifting_classification, scales = "free_y")

dev.off()

pdf("Shifting_ClassCount.pdf", width = 6, height = 8)

cdat_shifting %>%
  left_join(cdat_shifting_metrics) %>%
  filter(!duplicated(Stream_Name)) %>%
  group_by(shifting_classification) %>%
  tally() %>%
  mutate(shifting_classification = factor(shifting_classification, 
                                          levels=c("stable/mostly stable", "progressive",
                                                   "cyclical","non-specific cyclical",
                                                   "non-patterned change"))) %>%
  ggplot(aes(shifting_classification, n))+geom_bar(stat = "identity", fill="black")+
  theme_classic()+theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="", y="Count")

dev.off()

