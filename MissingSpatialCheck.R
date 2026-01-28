setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/SiSyn/ESOM")

sites<-read.csv("ESOM_Input_Sites_location.csv")

ref_table<-read.csv("Site_Reference_Table_01282026.csv", na.strings = "")

sites<-left_join(sites, ref_table[,c(3,23)])

sites<-sites %>%
  filter(!duplicated(Stream_Name))

sum(is.na(sites$Latitude))
sum(is.na(sites$Shapefile_Name))

sites_missing<-sites %>%
  filter(is.na(Shapefile_Name))

write.csv(sites_missing, "ESOM_MissingSpatialData.csv")

missing_latlong<-sites %>%
  filter(is.na(Latitude))

write.csv(missing_latlong, "Missing_LatLong.csv")
