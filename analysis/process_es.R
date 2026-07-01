#process
library(sf)
library(DT)
library(dplyr)
library(ggplot2)
library(terra)
library(SSDM)
library(purrr)

source("analysis/fkt_utils.R")

in_dir<-"data"
out_dir<-"output"

#study areas
studID<-c("SK021","FRL04","TRD")
#studID<-c("TRD")


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
    dplyr::left_join(es_mapping %>% dplyr::select(esID, userID, imp_acc, confidence),
                     by = c("esID", "userID"))




  esids<-unique(ind_pols$esID)

  ## site specific predictors
  pred<-SSDM::load_var(path=paste0(in_dir,"/",studID[i],"/2_env_var"), categorical = "lulc")


  algo_comp<-run_algorithm_comparison(es_ids = esids,
                                      studyID = studID[i],
                                      ind_pols = ind_pols,
                                      pred = pred,
                                      A_roi = A_roi,
                                      all_back_pts = all_back_pts,
                                      out_dir = paste0(out_dir,"/alg_performance"))


  run_participant_holdout(es_ids = esids,
                                      studyID = studID[i],
                                      ind_pols = ind_pols,
                                      pred = pred,
                                      A_roi = A_roi,
                                      all_back_pts = all_back_pts,
                                      out_dir = paste0(out_dir,"/participants_results"))



}

calc_algorithm_var(
  in_dir = "output/alg_performance"
)

calc_participant_var(
  in_dir = "output/participants_results",
  out_dir=paste0(in_dir, "/part_var")
)



# explain participant / algorithm uncertainty
df_part <- build_uncertainty_df (
  part_uncertainty_dir = "output/participants_results/part_var",
  alg_dir = "output/alg_performance",
  env_base_dir = "data")

#
df_part$studyID <- factor(df_part$studyID)
df_part$ES <- factor(df_part$ES)
df_part$lulc<- trunc(df_part$lulc / 100)
df_part<-df_part%>%filter(lulc != 0)
df_part$lulc <- factor(df_part$lulc)

saveRDS(df_part,"output/processed/df_part.rds")


library(mgcv)


## model per ES
es_ids <- sort(unique(df_part$ES))

gam_summary <- vector("list", length(es_ids))

for(i in seq_along(es_ids)){

  es <- es_ids[i]

  message("Processing ", es)

  dat <- df_part %>%
    filter(ES == es)

  m <- gam(
    R ~
      s(dem) +
      s(acc) +
      s(int) +
      factor(lulc) +
      factor(studyID),
    data = dat,
    family = betar(),
    method = "REML"
  )

  sm <- summary(m)

  ## Smooth statistics
  s_tab <- as.data.frame(sm$s.table)

  p_tab <- as.data.frame(sm$p.table)

  ## helper function
  coef_or_na <- function(x){
    if(x %in% rownames(p_tab))
      p_tab[x,"Estimate"]
    else
      NA
  }

  gam_summary[[i]] <- data.frame(
    ES = es,
    n = nrow(dat),
    Adj_R2 = sm$r.sq,
    expl_dev = sm$dev.expl,
    REML = sm$sp.criterion,
    F_dem = s_tab["s(dem)", "F"],
    edf_dem = s_tab["s(dem)", "edf"],
    F_acc = s_tab["s(acc)", "F"],
    edf_acc = s_tab["s(acc)", "edf"],
    F_int = s_tab["s(int)", "F"],
    edf_int = s_tab["s(int)", "edf"],
    Delta_SK021 = coef_or_na("studyIDSK021"),
    Delta_TRD   = coef_or_na("studyIDTRD")
  )
}

gam_summary <- bind_rows(gam_summary)

write.csv(
  gam_summary,
  "output/GAM_ES_summary_R.csv",
  row.names = FALSE
)


  es_ids <- sort(unique(df_part$ES))

  out <- list()
  for(es in es_ids){

    message(es)

    dat <- filter(df_part, ES == es)
    dat <- dat %>%
      filter(lulc != 0) %>%
      droplevels()


    m <- gam(
      R ~
        s(dem, by = studyID) +
        s(acc, by = studyID) +
        s(int, by = studyID) +
        lulc +
        studyID,
      data = dat,
      family = betar(),
      method = "REML"
    )

    message("model done")


    pd <- dat %>%
      split(.$studyID) %>%
      imap_dfr(function(d, id) {

        bind_rows(
          calc_pdp(d, m, "dem") %>%
            mutate(predictor = "Elevation"),

          calc_pdp(d, m, "acc") %>%
            mutate(predictor = "Accessibility"),
          #
          calc_pdp(d, m, "int") %>%
            mutate(predictor = "Integrity")
        ) %>%
          mutate(studyID = id)
      })
    out[[es]] <- pd %>%
      mutate(ES = es)

    ggplot(pd,
           aes(x, R,
               colour = studyID,
               fill = studyID)) +
      geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        colour = "grey50",
        linewidth = 0.5
      ) +
      geom_ribbon(aes(ymin = pmax(0, R - SD),
                      ymax = pmin(1, R + SD)),
                  alpha = 0.15,
                  colour = NA) +
      geom_line(linewidth = 1) +
      facet_wrap(~predictor, scales = "free_x") +
      theme_classic()

}

  pd_all <- bind_rows(out)
  trd<-pd_all%>%filter(studyID == "TRD")

  df <- pd_all %>%
    mutate(
      x = case_when(
        predictor == "Integrity" & studyID == "TRD" ~ x / 100,
        TRUE ~ x
      )
    )

  saveRDS(
    pd_all,
    "output/pdp_all_ES2.rds"
  )
