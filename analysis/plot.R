# plots

library(sf)
library(DT)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(terra)
source("analysis/fkt_utils.R")


## plot map


## plot algo performance per ES per studArea
## Read all algorithm performance tables
files <- list.files(
  "output/alg_performance",
  pattern = "_algorithm_performance\\.rds$",
  full.names = TRUE
)

alg_perf <- lapply(files, readRDS) |>
  bind_rows()
names(alg_perf)

ggplot(alg_perf,
       aes(x = Algorithm,
           y = ES,
           fill = AUC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", AUC)),
            size = 3) +
  facet_wrap(~studyID) +
  scale_fill_viridis_c(
    name = "AUC",
    limits = c(0.5, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    panel.grid = element_blank()
  )


### plot 2 var imp
files <- list.files(
  "output/alg_performance",
  pattern = "_importance\\.rds$",
  full.names = TRUE
)

imp <- bind_rows(lapply(files, readRDS))

## Keep only RF
imp_rf <- imp %>%
  filter(Algorithm == "RF")
imp_rank <- imp_rf %>%
  group_by(studyID, ES, Algorithm) %>%
  mutate(
    Rank = rank(-Importance, ties.method = "first")
  ) %>%
  ungroup()




ggplot(imp_rank,
       aes(x = Variable,
           y = ES,
           fill = Rank)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(Rank,0)), size = 4) +
  facet_wrap(~studyID) +
  scale_fill_gradient(
    low = "darkgreen",
    high = "white",
    trans = "reverse",
    breaks = 1:4
  ) +
  labs(
    fill = "Rank",
    x = "Predictor",
    y = "Algorithm"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


## participant uncertainty
# part_unc<-terra::rast("output/participants_results/FRL04_aest_userStack.tif")
# plot(part_unc)
#Mean absolute deviation from the full model after removing one participant.


files <- list.files(
  "output/participants_results",
  pattern = "_participant_performance\\.rds$",
  full.names = TRUE
)

part_unc <- bind_rows(lapply(files, readRDS))

### test clustering:
# test<-part_unc%>%filter(studyID == "FRL04",ES == "aest")
# test <- test |>
#   dplyr::select(
#     delta_AUC,
#     MeanDiff,
#     RMSE,
#     Corr,
#     Delta_dem,
#     Delta_lulc,
#     Delta_int,
#     Delta_acc
#   )
#
# X <- scale(test)
#
#
# hc <- hclust(dist(X), method = "ward.D2")
# plot(hc)
#
# clusters <- cutree(hc, k = 3)
#
# test$cluster <- factor(clusters)



ggplot(part_unc,
       aes(x = ES,
           y = MeanDiff)) +
  geom_boxplot(fill = "grey80",
               outlier.shape = NA) +
  geom_jitter(width = 0.15,
              alpha = 0.7,
              size = 2) +
  facet_wrap(~studyID, nrow = 1) +
  labs(
    x = "Ecosystem service",
    y = "Mean leave-one-participant prediction difference"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


##  R 0=algorithm uncertainty higher, R = 1 participant uncertainty higher

mean_var_ratio_plot<-plot_variance_ratio("output/alg_performance/var_alg",
                                         "output/participants_results")
plot(mean_var_ratio_plot$FRL04)



## relationship algo and mapper (standardized)
# part_sum<-part_unc%>%group_by(studyID,ES)%>%summarise(MeanPartDiff = mean(MeanDiff))
# part_sum<-part_sum%>%select(studyID,ES,MeanPartDiff)
# colnames(part_sum)<-c("studyID","ES","u_participant")

part_var_pix<-read_var_summary(in_dir="output/participants_results/part_var/",suffix = "_V_part")
part_var_pix$ES <- gsub("_V_part", "", part_var_pix$ES)
colnames(part_var_pix)<-c("studyID","ES","var_participant")

# algo uncertainty from raster (also standardized):
algo_var_pix<-read_uncertainty_summary("output/alg_performance/var_alg",suffix = "_V_alg")
algo_var_pix$ES <- gsub("_V_alg", "", algo_var_pix$ES)
colnames(algo_var_pix)<-c("studyID","ES","var_algorithm")

# join with studyID and ES
uncertainty_all<-merge(part_var_pix,algo_var_pix,by=c("studyID","ES"))

uncertainty_all <- uncertainty_all %>%
  dplyr::group_by(studyID) %>%
  dplyr::mutate(
    Alg = var_algorithm / max(var_algorithm),
    Part = var_participant / max(var_participant)
  ) %>%
  dplyr::ungroup()

##


ggplot(
  uncertainty_all,
  aes(
    x = var_algorithm,
    y = var_participant,
    colour = ES
  )
) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "grey70",
    linewidth = 0.6
  ) +


  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +

  ggrepel::geom_text_repel(
    aes(label = ES),
    size = 3,
    max.overlaps = 6,
    show.legend = FALSE
  ) +

  facet_wrap(~studyID) +

  coord_equal() +

  scale_colour_viridis_d(option = "D") +

  labs(
    x = "Algorithm variation",
    y = "Participant variation",
    colour = "Ecosystem service"
  ) +

  theme_minimal(base_size = 13) +

  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    aspect.ratio = 1
  )



ggplot(df,
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
  #facet_wrap(~predictor+ES, scales = "free_x") +
  facet_grid(vars(ES), vars(predictor),scales = "free_x" )+
  theme_classic()
