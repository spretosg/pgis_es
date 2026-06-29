#process
library(sf)
library(DT)
library(dplyr)
library(ggplot2)
library(terra)
library(SSDM)

source("analysis/fkt_utils.R")

in_dir<-"data"
out_dir<-"output"

#study areas
# studID<-c("SK021","FRL04","TRD")
studID<-c("TRD")


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
  ind_pols<-read_sf(paste0(in_dir,"/",studID[i],"/ind_polys_R1.gpkg"))%>%dplyr::filter(siteID == "rescape_TRD")

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


  run_participant_holdout(es_ids = esids,
                                      studyID = studID[i],
                                      mappers = users,
                                      ind_pols = ind_pols,
                                      pred = pred,
                                      A_roi = A_roi,
                                      all_back_pts = all_back_pts,
                                      out_dir = paste0(out_dir,"/participants_results"))



}

calc_algorithm_uncertainty(
  "output",
  file.path(in_dir, "algo_uncertainty")
)

calc_participant_uncertainty(
  in_dir = "output/participants_results",
  file.path(in_dir, "part_uncertainty")
)


library(terra)
library(ggplot2)

plot_uncertainty_ratio <- function(
    alg_dir,
    part_dir){

  alg_files <- list.files(
    alg_dir,
    pattern = "_MeanAlgorithmUncertainty\\.tif$",
    full.names = TRUE
  )

  plot_list <- list()

  for(f_alg in alg_files){

    studyID <- sub(
      "_MeanAlgorithmUncertainty\\.tif",
      "",
      basename(f_alg)
    )

    f_part <- file.path(
      part_dir,
      paste0(studyID,
             "_MeanParticipantUncertainty.tif")
    )

    if(!file.exists(f_part)){
      warning(studyID, ": participant map not found")
      next
    }

    Ualg  <- rast(f_alg)
    Upart <- rast(f_part)

    ## ratio
    R <- Upart / (Upart + Ualg)

    names(R) <- "Ratio"

    ## convert for ggplot
    d <- as.data.frame(R,
                       xy = TRUE,
                       na.rm = TRUE)

    p <- ggplot(d,
                aes(x, y, fill = Ratio)) +
      geom_raster() +
      coord_equal() +
      scale_fill_gradient2(
        low = "#2b83ba",
        mid = "white",
        high = "#d7191c",
        midpoint = 0.5,
        limits = c(0,1),
        name = expression(R)
      ) +
      labs(title = studyID) +
      theme_void() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )

    plot_list[[studyID]] <- p

  }

  return(plot_list)

}

