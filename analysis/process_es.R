#process
library(sf)
library(leaflet)
library(DT)
library(dplyr)
library(ggplot2)
library(dplyr)
library(terra)
library(SSDM)

source("analysis/fkt_utils.R")

in_dir<-"data"
out_dir<-"output"

#study areas
studID<-c("SK021","FRL04","TRD")

#load stud areas
stud_site<-read_sf(paste0(in_dir,"/studArea.gpkg"))



for(i in 1:length(studID)){
  #
  tmpSite<-stud_site%>%filter(siteID == studID[i])
  A_roi<-tmpSite$siteAREAkm2*10^6
  resolution = 250^2
  all_back_pts<- round(A_roi/resolution,0)

  ## join imp_access to ind_pols via esID and userID from esmappingR1.csv
  es_mapping<-read.csv(paste0(in_dir,"/",studID[i],"/es_mappingR1.csv"))
  ind_pols<-read_sf(paste0(in_dir,"/",studID[i],"/ind_polys_R1.gpkg"))%>%dplyr::filter(siteID == studID[i])

  ind_pols <- ind_pols %>%
    dplyr::left_join(es_mapping %>% dplyr::select(esID, userID, imp_acc),
                     by = c("esID", "userID"))

  users<-unique(ind_pols$userID)
  esids<-unique(ind_pols$esID)

  ## site specific predictors
  pred<-SSDM::load_var(path=paste0(in_dir,"/",studID[i],"/2_env_var"), categorical = "lulc")


  algo_comp<-run_algorithm_comparison(es_ids = esids,
                                      studyID = studID[i],
                                      ind_pols = ind_pols,
                                      pred = pred,
                                      A_roi = A_roi,
                                      all_back_pts = all_back_pts,
                                      out_dir = out_dir)

}

