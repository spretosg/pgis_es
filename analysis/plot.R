# plots

library(sf)
library(DT)
library(dplyr)
library(ggplot2)
library(terra)
source("analysis/fkt_utils.R")


## plot map


## plot algo performance per ES per studArea
## Read all algorithm performance tables
files <- list.files(
  "output",
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


### plot 2 var imp for RF
files <- list.files(
  "output",
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
part_sum<-part_unc%>%group_by(studyID,ES)

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

plots <- plot_uncertainty_ratio(
  alg_dir = "output/uncertainty",
  part_dir = "output/part_uncertainty"
)

plot(plots$SK021)

## relationship algo and mapper
part_sum<-part_unc%>%group_by(studyID,ES)%>%summarise(MeanPartDiff = mean(MeanDiff))
part_sum<-part_sum%>%select(studyID,ES,MeanPartDiff)
colnames(part_sum)<-c("studyID","ES","u_participant")

# algo uncertainty from raster:
algo_unc_pix<-read_uncertainty_summary("output/uncertainty/","_algo_Ualgorithm")
algo_unc_pix$ES <- gsub("_algo_Ualgorithm", "", algo_unc_pix$ES)
colnames(algo_unc_pix)<-c("studyID","ES","u_algorithm")

# join with studyID and ES
uncertainty_all<-merge(part_sum,algo_unc_pix,by=c("studyID","ES"))

##

ggplot(uncertainty_all,
       aes(x = u_algorithm,
           y = u_participant,
           label = ES)) +

  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2,
              colour = "grey50") +

  geom_point(aes(colour = ES),
             size = 3) +

  ggrepel::geom_text_repel(size = 3,
                           show.legend = FALSE) +

  facet_wrap(~studyID) +

  coord_equal() +

  labs(
    x = "Mean algorithm uncertainty",
    y = "Mean participant uncertainty"
  ) +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
