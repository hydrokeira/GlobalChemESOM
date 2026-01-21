require(data.table)
require(readxl)
require(dplyr)
require(tidyr)
require(stringr)
require(ggplot2)

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

ref_table<-read_xlsx("Site_Reference_Table_01092026.xlsx")

good_sites_latlong<-left_join(good_sites, ref_table[,c(1,3,9,10)])

good_sites_latlong<-good_sites_latlong[!duplicated(good_sites_latlong$Stream_Name),]

good_sites_latlong <- good_sites_latlong %>%
  filter(!Stream_Name=="Sopchoppy River")

unique(good_sites_latlong$LTER)

write.csv(good_sites_latlong, "ESOM_Input_Sites_location.csv")

chem_cast_crop_CC <- chem_cast_crop_CC %>%
  filter(Stream_Name %in% good_sites$Stream_Name)

chem_byYear<-chem_cast_crop_CC %>%
  mutate(date=as.Date(date)) %>%
  group_by(year(date)) %>%
  tally()

chem_byYear %>%
  filter(!is.na(n)) %>%
  ggplot(aes(x=`year(date)`, y=1, fill=n, color=n))+geom_tile()+
  scale_fill_gradient(low="grey", high = "black")+
  scale_color_gradient(low="grey", high = "black")+
  theme_classic()+labs(x="Year", y="",fill="Number of observations",
                       color="Number of observations")+
  theme(text = element_text(size=20), legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

pdf("Num_Obs_Year.pdf", width = 10, height = 3)

chem_byYear %>%
  filter(!is.na(n)) %>%
  ggplot(aes(x=`year(date)`, y=n))+geom_bar(stat="identity", col="white", fill="black")+
  #scale_fill_gradient(low="grey", high = "black")+
  #scale_color_gradient(low="grey", high = "black")+
  theme_classic()+labs(x="Year",y="Number of observations")+
  theme(text = element_text(size=15), legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1))

dev.off()

chem_byYear_bySite<-chem_cast_crop_CC %>%
  mutate(date=as.Date(date), year=year(date)) %>%
  dplyr::group_by(Stream_Name) %>%
  summarise(n_years=n_distinct(year))

pdf("Num_Years_Site.pdf", width = 10, height = 3)

chem_byYear_bySite %>%
  ggplot(aes(x=n_years))+geom_bar(stat="count", col="white", fill="black")+
  #scale_fill_gradient(low="grey", high = "black")+
  #scale_color_gradient(low="grey", high = "black")+
  theme_classic()+labs(x="Number of Years in Period of Record",y="Number of Sites")+
  theme(text = element_text(size=15), legend.position = "bottom", 
        legend.text = element_text(angle = 45, hjust = 1))

dev.off()
