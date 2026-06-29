#functions


## function to sample presence points
make_presence_points <- function(polys, A_roi, all_back_pts, min_in_pts = 10) {

  pts_list <- lapply(1:nrow(polys), function(i) {

    A_tmp <- as.numeric(st_area(polys[i,]))
    prop <- A_tmp / A_roi
    npts <- round(all_back_pts * prop)
    npts <- max(npts, min_in_pts)

    # ES value scaling
    es_scale <- (1 + polys$es_value[i]) / 5
    npts <- round(npts * es_scale)

    # sample points
    pts_sf <- st_sample(polys[i,], npts, type = "random") |>
      st_as_sf() |>
      mutate(inside = 1)

    pts_sf
  })

  do.call(rbind, pts_list)
}

## algorithm performance for studyID, esID
run_algorithm_comparison <- function(es_ids,
                                     studyID,
                                     ind_pols,
                                     pred,
                                     A_roi,
                                     all_back_pts,
                                      out_dir = "algorithm_results",
                                     algos = c("GLM","RF","GAM","MAXENT")) {


  performance_all <- list()
  imp_all<-list()

  for (es in es_ids) {

    message("Processing ES ", es)

    pol_es <- dplyr::filter(ind_pols, esID == es)

    if (nrow(pol_es) == 0)
      next
    pol_es<-st_make_valid(pol_es)
    ## Presence points
    pts_full <- make_presence_points(
      pol_es,
      A_roi,
      all_back_pts,
      min_in_pts = 10
    )

    pts_full <- sf::st_transform(pts_full, sf::st_crs(pred))
    pts_full_sp <- as(pts_full, "Spatial")

    ## Predictor stack
    pred_w <- raster::stack(
      pred$dem * 1,
      pred$lulc * 1,
      pred$int * 1,
      pred$acc * (mean(pol_es$imp_acc, na.rm = TRUE) / 5)
    )

    extracted <- raster::extract(pred_w, pts_full_sp)

    df_full <- cbind(pts_full, extracted)
    df_full <- cbind(df_full, sf::st_coordinates(df_full))

    colnames(df_full)[colnames(df_full) %in% c("X","Y")] <- c("lon","lat")

    df_full <- sf::st_drop_geometry(df_full)
    df_full <- na.omit(df_full)
    df_full <- dplyr::select(df_full, lon, lat)

    df_full$SPECIES <- "pres"

    ## Fit models
    mods <- lapply(algos, function(a)
      SSDM::modelling(
        a,
        df_full,
        pred_w,
        Xcol = "lon",
        Ycol = "lat",
        cv = "holdout",
        cv.param = c(0.7,2),
        final.fit.data = "all"
      )
    )

    names(mods) <- algos

    ## Performance table
    perf <- do.call(rbind,
                    lapply(names(mods), function(a){

                      p <- as.data.frame(mods[[a]]@evaluation)

                      p$Algorithm <- a
                      p$studyID <- studyID
                      p$ES <- es

                      p

                    })
    )

    performance_all[[as.character(es)]] <- perf

    ## Prediction stack
    alg_stack <- terra::rast(
      lapply(mods, function(m)
        terra::rast(m@projection))
    )

    names(alg_stack) <- names(mods)

    terra::writeRaster(
      alg_stack,
      filename = file.path(
        out_dir,
        paste0(studyID, "_", es, "_algo.tif")
      ),
      overwrite = TRUE
    )

    ## Variable importance
    imp <- do.call(rbind,
                   lapply(names(mods), function(a){

                     x <- mods[[a]]@variable.importance
                     x <- x/sum(x)

                     data.frame(
                       studyID = studyID,
                       ES = es,
                       Algorithm = a,
                       Variable = names(x),
                       Importance = as.numeric(x)
                     )

                   })
    )
    imp_all[[as.character(es)]] <- imp



  }

  imp_all <- dplyr::bind_rows(imp_all)

  performance_all <- dplyr::bind_rows(performance_all)
  row.names(performance_all)<-NULL

  saveRDS(
    imp_all,
    file.path(
      out_dir,
      paste0(studyID, "_importance.rds")
    )
  )

  saveRDS(
    performance_all,
    file.path(
      out_dir,
      paste0(studyID, "_algorithm_performance.rds")
    )
  )

  return(performance_all)

}

run_participant_holdout <- function(es_ids,
                                     studyID,
                                     mappers,
                                     ind_pols,
                                     pred,
                                     A_roi,
                                     all_back_pts,
                                     out_dir = "participants_results") {

  performance_all <- list()
  imp_all<-list()

  for (es in es_ids) {

    message("Processing ES ", es)

    pol_es <- dplyr::filter(ind_pols, esID == es)

    if (nrow(pol_es) == 0)
      next
    pol_es<-st_make_valid(pol_es)

    ## Presence points
    pts_full <- make_presence_points(
      pol_es,
      A_roi,
      all_back_pts,
      min_in_pts = 10
    )

    pts_full <- sf::st_transform(pts_full, sf::st_crs(pred))
    pts_full_sp <- as(pts_full, "Spatial")

    ## Predictor stack
    pred_w <- raster::stack(
      pred$dem * 1,
      pred$lulc * 1,
      pred$int * 1,
      pred$acc * (mean(pol_es$imp_acc, na.rm = TRUE) / 5)
    )

    extracted <- raster::extract(pred_w, pts_full_sp)

    df_full <- cbind(pts_full, extracted)
    df_full <- cbind(df_full, sf::st_coordinates(df_full))

    colnames(df_full)[colnames(df_full) %in% c("X","Y")] <- c("lon","lat")

    df_full <- sf::st_drop_geometry(df_full)
    df_full <- na.omit(df_full)
    df_full <- dplyr::select(df_full, lon, lat)

    df_full$SPECIES <- "pres"

    m_full <- SSDM::modelling(
      "RF",
      df_full,
      pred_w,
      Xcol = "lon",
      Ycol = "lat",
      cv = "holdout",
      cv.param = c(0.7,2),
      final.fit.data = "all"
    )

    auc_full <- m_full@evaluation$AUC
    varimp_full <- m_full@variable.importance
    proj_full <- terra::rast(m_full@projection)






    ## Table to store results ------------------------------------------------

    user_scores <- data.frame(
      studyID = studyID,
      ES = es,
      userID = mappers,
      AUC_full = auc_full,
      AUC_minusUser = NA_real_,
      delta_AUC = NA_real_,
      MeanDiff = NA_real_,
      RMSE = NA_real_,
      Corr = NA_real_,
      Delta_dem = NA_real_,
      Delta_lulc = NA_real_,
      Delta_int = NA_real_,
      Delta_acc = NA_real_
    )

    ## Store prediction maps if desired
    projection_list <- list()
    projection_results<-list()
    importance_all <- list()

    ## Leave-one-user-out ----------------------------------------------------

    for (u in mappers) {

      message("Removing user: ", u)

      pol_minus <- pol_es %>%
        filter(userID != u)

      pts_minus <- make_presence_points(pol_minus, A_roi, all_back_pts)
      pts_minus <- st_transform(pts_minus, st_crs(pred))
      pts_minus_sp <- as(pts_minus, "Spatial")

      pred_w2 <- raster::stack(
        pred$dem * 1,
        pred$lulc * 1,
        pred$int * 1,
        pred$acc * (mean(pol_minus$imp_acc, na.rm = TRUE) / 5)
      )

      extracted2 <- raster::extract(pred_w2, pts_minus_sp)

      df_minus <- cbind(pts_minus, extracted2)
      df_minus <- cbind(df_minus, st_coordinates(df_minus))
      colnames(df_minus)[colnames(df_minus) %in% c("X","Y")] <- c("lon","lat")

      df_minus <- st_drop_geometry(df_minus)
      df_minus <- na.omit(df_minus)
      df_minus <- dplyr::select(df_minus, lon, lat)

      df_minus$SPECIES <- "pres"

      ## Fit model
      m_minus <- SSDM::modelling(
        "RF",
        df_minus,
        pred_w,
        Xcol="lon",
        Ycol="lat",
        cv="holdout",
        cv.param=c(0.7,2),
        final.fit.data="all"
      )

      ## Evaluation
      auc_minus <- m_minus@evaluation$AUC

      ## Variable importance
      varimp_minus <- m_minus@variable.importance

      ## Prediction map
      proj_minus <- terra::rast(m_minus@projection)

      projection_list[[as.character(u)]] <- proj_minus

      ## Spatial difference
      diff <- proj_full - proj_minus

      mean_diff <- terra::global(abs(diff), "mean", na.rm=TRUE)[1,1]

      rmse <- sqrt(
        terra::global(diff^2, "mean", na.rm=TRUE)[1,1]
      )

      vals <- na.omit(cbind(
        terra::values(proj_full),
        terra::values(proj_minus)
      ))

      corr <- cor(vals[,1],
                  vals[,2],
                  use = "complete.obs")

      ## Save

      user_scores[user_scores$userID==u,"AUC_minusUser"] <- auc_minus
      user_scores[user_scores$userID==u,"delta_AUC"] <- auc_full-auc_minus

      user_scores[user_scores$userID==u,"MeanDiff"] <- mean_diff
      user_scores[user_scores$userID==u,"RMSE"] <- rmse
      user_scores[user_scores$userID==u,"Corr"] <- corr

      user_scores[user_scores$userID==u,"Delta_dem"] <-
        varimp_full[1]-varimp_minus[1]

      user_scores[user_scores$userID==u,"Delta_lulc"] <-
        varimp_full[2]-varimp_minus[2]

      user_scores[user_scores$userID==u,"Delta_int"] <-
        varimp_full[3]-varimp_minus[3]

      user_scores[user_scores$userID==u,"Delta_acc"] <-
        varimp_full[4]-varimp_minus[4]

      importance_all[[u]] <- data.frame(
        studyID = studyID,
        ES = es,
        userID = u,
        dem = varimp_minus[1],
        lulc = varimp_minus[2],
        int = varimp_minus[3],
        acc = varimp_minus[4]
      )





    }#/user loop

    # imp_es <- data.frame(
    #   studyID = studyID,
    #   ES = es,
    #   userID = users_es,
    #   dem = user_scores$Delta_dem,
    #   lulc = user_scores$Delta_lulc,
    #   int = user_scores$Delta_int,
    #   acc = user_scores$Delta_acc
    # )

    imp_es <- dplyr::bind_rows(importance_all)
    performance_all[[as.character(es)]] <- user_scores
    imp_all[[as.character(es)]] <- imp_es
    user_stack <- terra::rast(projection_list)

    names(user_stack) <- names(projection_list)

    terra::writeRaster(
      user_stack,
      filename = file.path(
        out_dir,
        paste0(studyID,"_", es, "_userStack.tif")
      ),
      overwrite = TRUE
    )

    # results_list[[as.character(es)]] <- user_scores
    # projection_results[[as.character(es)]] <- projection_list

        # preds <- rast(projection_results[[as.character(es)]])
    #
    # names(preds) <- names(projection_results[[as.character(es)]])


  } #/es loop
    performance_all <- dplyr::bind_rows(performance_all)
    imp_all <- dplyr::bind_rows(imp_all)

  saveRDS(
    performance_all,
    file.path(
      out_dir,
      paste0(studyID, "_participant_performance.rds")
    )
  )

  saveRDS(
    imp_all,
    file.path(
      out_dir,
      paste0(studyID, "_participant_importance.rds")
    )
  )

  return(list(
    performance = performance_all,
    importance = imp_all
  ))


}

library(terra)
library(stringr)

calc_algorithm_uncertainty <- function(
    in_dir,
    out_dir
){

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(
    in_dir,
    pattern = "_algo\\.tif$",
    full.names = TRUE
  )

  study_list <- list()

  for(f in files){

    message("Processing ", basename(f))

    nm <- tools::file_path_sans_ext(basename(f))
    nm <- sub("_algo$", "", nm)

    studyID <- strsplit(nm, "_")[[1]][1]
    ES <- paste(strsplit(nm, "_")[[1]][-1], collapse="_")

    stk <- rast(f)

    ## pixel-wise uncertainty
    Ualg <- app(stk, sd, na.rm = TRUE)

    names(Ualg) <- "U_algorithm"

    writeRaster(
      Ualg,
      file.path(
        out_dir,
        paste0(studyID, "_", ES, "_Ualgorithm.tif")
      ),
      overwrite=TRUE
    )

    if(is.null(study_list[[studyID]]))
      study_list[[studyID]] <- list()

    study_list[[studyID]][[ES]] <- Ualg
  }

  ## average across ES for each study
  for(st in names(study_list)){

    es_stack <- rast(study_list[[st]])

    Umean <- app(es_stack, mean, na.rm=TRUE)

    names(Umean) <- "Mean_Ualgorithm"

    writeRaster(
      Umean,
      file.path(
        out_dir,
        paste0(st, "_MeanAlgorithmUncertainty.tif")
      ),
      overwrite=TRUE
    )

  }

  invisible(study_list)

}

library(terra)

calc_participant_uncertainty <- function(
    in_dir,
    out_dir
){

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(
    in_dir,
    pattern = "_userStack\\.tif$",
    full.names = TRUE
  )

  study_list <- list()

  for(f in files){

    message("Processing ", basename(f))

    nm <- tools::file_path_sans_ext(basename(f))
    nm <- sub("_userStack$", "", nm)

    parts <- strsplit(nm, "_")[[1]]

    studyID <- parts[1]
    ES <- paste(parts[-1], collapse = "_")

    stk <- rast(f)

    Upart <- app(stk, sd, na.rm = TRUE)

    names(Upart) <- "U_participant"

    writeRaster(
      Upart,
      file.path(
        out_dir,
        paste0(studyID, "_", ES, "_Uparticipant.tif")
      ),
      overwrite = TRUE
    )

    if(is.null(study_list[[studyID]]))
      study_list[[studyID]] <- list()

    study_list[[studyID]][[ES]] <- Upart
  }

  ## Mean uncertainty across ES
  for(st in names(study_list)){

    es_stack <- rast(study_list[[st]])

    Umean <- app(es_stack, mean, na.rm = TRUE)

    names(Umean) <- "Mean_Uparticipant"

    writeRaster(
      Umean,
      file.path(
        out_dir,
        paste0(st, "_MeanParticipantUncertainty.tif")
      ),
      overwrite = TRUE
    )

  }

  invisible(study_list)

}

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


library(terra)
library(dplyr)

read_uncertainty_summary <- function(in_dir,
                                     suffix){

  files <- list.files(
    in_dir,
    pattern = paste0(suffix, "\\.tif$"),
    full.names = TRUE
  )

  res <- lapply(files, function(f){

    r <- rast(f)

    nm <- tools::file_path_sans_ext(basename(f))
    nm <- sub(paste0("_", suffix, "$"), "", nm)

    parts <- strsplit(nm, "_")[[1]]

    studyID <- parts[1]
    ES <- paste(parts[-1], collapse = "_")

    data.frame(
      studyID = studyID,
      ES = ES,
      MeanUncertainty = terra::global(r, "mean", na.rm = TRUE)[1,1]
    )

  })

  bind_rows(res)

}
