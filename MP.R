### ------------------------------------------------------------------------ ###
### some seasonal projections ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### set-up ####
### ------------------------------------------------------------------------ ###
source("funs.R")
library(FLasher)
library(doParallel)
library(tidyverse)
library(cowplot)

cl <- makeCluster(10)
registerDoParallel(cl)
. <- foreach(seq(length(cl))) %dopar% {
  source("funs.R")
  library(FLasher)
}

### ------------------------------------------------------------------------ ###
### HR with annual TAC - lags ####
### ------------------------------------------------------------------------ ###

stk_rnd_500 <- readRDS("input/san/stk_rnd_500.rds")
sr_500 <- readRDS("input/san/sr_500.rds")
hr_res_MSY_annual <- readRDS("input/san/hr_res_MSY_annual.rds")
hr_target <- hr_res_MSY_annual$maximum ### annual MSY HR
refpts <- readRDS("input/san/refpts.rds")

iters_list <- split(seq(500), sort(seq(500) %% 10))
input_list <- lapply(iters_list, function(x) {
  list(stk = iter(stk_rnd_500, x), sr = iter(sr_500, x))
})

### run for lag 0, 1, ... 8
stk_lags <- foreach(lag = 0:8) %do% {
  print(paste0("lag=", lag))
  stk_junks <- foreach(input = input_list) %dopar% {
    # browser()
    res_i <- mse_loop(MP = "hr", force_seasonal = TRUE,
                      om_stk = input$stk, om_sr = input$sr, 
                      yrs = 101:200, seasons = 1:4, 
                      target = hr_target, lag = lag, catch_interval = 4,
                      verbose = FALSE)
    return(res_i)
  }
  stk_out <- stk_rnd_500
  for (i in seq_along(stk_junks)) 
    iter(stk_out, iters_list[[i]]) <- stk_junks[[i]]
  saveRDS(stk_out, paste0("output/san/seasonal/rnd_MSY_hr_annual_lag", lag,
                          ".rds"))
}
plot(stk_out)
plot(window(stk_out, start = 121, end = 125), iter = 1:5)

### create list with all results
stk_lags <- lapply(0:8, function(lag) {
  readRDS(paste0("output/san/seasonal/rnd_MSY_hr_annual_lag", lag, ".rds"))
})
names(stk_lags) <- 0:8

stats <- lapply(seq_along(stk_lags), function(x) {
  #browser()
  qnts <- FLQuants("Catch" = seasonSums(catch(stk_lags[[x]])[, ac(101:125)])/
                     refpts["msy", "yield"],
                   "Fbar" = seasonSums(fbar(stk_lags[[x]])[, ac(101:125)])/
                     refpts["msy", "harvest"],
                   "SSB" = ssb(stk_lags[[x]])[, ac(101:125),, 1]/
                     refpts["msy", "ssb"],
                   "TSB" = tsb(stk_lags[[x]])[, ac(101:125),, 1]/
                     refpts["virgin", "biomass"],
                   "risk" = ssb(stk_lags[[x]])[, ac(101:125),, 1] %=% NA_real_)
  qnts$risk[] <- as.numeric(ssb(stk_lags[[x]])[, ac(101:125),, 1] < 163)
  qnts <- as.data.frame(qnts)
  qnts$lag <- as.numeric(names(stk_lags)[x])
  return(qnts)
})
stats <- do.call(rbind, stats)

p_catch <- stats %>%
  filter(qname == "Catch") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(annual~catch/MSY), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 2.5)) +
  facet_wrap(~ "Harvest rate") + 
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_catch
p_ssb <- stats %>%
  filter(qname == "SSB") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(SSB/SSB[MSY]), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ssb
p_tsb <- stats %>%
  filter(qname == "TSB") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(B/B[0]), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_tsb
p_risk <- stats %>%
  filter(qname == "risk") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  #geom_hline(yintercept = 1, colour = "grey") +
  #geom_violin(fill = "grey", size = 0.2, show.legend = FALSE) +
  geom_col(data = stats %>%
               filter(qname == "risk") %>%
               group_by(lag) %>%
               summarise(data = mean(data)),
             aes(x = lag/4, y = data, group = lag), 
           fill = "grey", colour = "black", size = 0.2) +
  #scale_fill_brewer(palette = "Set1") +
  # geom_boxplot(fill = "white", width = 0.1, size = 0.2,
  #              outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
  #              outlier.fill = "transparent") +
  labs(y = expression(B[lim]~risk), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"))
p_risk

p_hr_stats <- plot_grid(p_catch, p_ssb, p_tsb, p_risk, ncol = 1, 
                        align = "v", rel_heights = c(1, 1, 1, 1.2))
p_hr_stats
ggsave("output/plots/san_MSY_HR_lag_stats.png",
       width = 8.5, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_HR_lag_stats.pdf",
       width = 8.5, height = 10, units = "cm")



### ------------------------------------------------------------------------ ###
### Escapement strategy with annual TAC - lags ####
### ------------------------------------------------------------------------ ###

stk_rnd_500 <- readRDS("input/san/stk_rnd_500.rds")
sr_500 <- readRDS("input/san/sr_500.rds")
esc_res_MSY_annual <- readRDS("input/san/esc_res_annual_MSY.rds")
esc_target <- esc_res_MSY_annual$maximum ### annual MSY HR
refpts <- readRDS("input/san/refpts.rds")

iters_list <- split(seq(500), sort(seq(500) %% 10))
input_list <- lapply(iters_list, function(x) {
  list(stk = iter(stk_rnd_500, x), sr = iter(sr_500, x))
})

### run for lag 0, 1, ... 8
stk_esc_lags <- foreach(lag = 0:8) %do% {
  print(paste0("starting lag=", lag))
  stk_junks <- foreach(input = input_list) %dopar% {
    # browser()
    res_i <- mse_loop(MP = "escapement", om_stk = input$stk, om_sr = input$sr, 
                      yrs = 101:200, seasons = 1:4, 
                      target = esc_target, lag = lag, catch_interval = 4,
                      verbose = FALSE)
    return(res_i)
  }
  stk_out <- stk_rnd_500
  for (i in seq_along(stk_junks)) 
    iter(stk_out, iters_list[[i]]) <- stk_junks[[i]]
  saveRDS(stk_out, paste0("output/san/seasonal/rnd_MSY_esc_annual_lag", lag,
                          ".rds"))
  print(paste0("finished lag=", lag))
}
### create list with all results
stk_esc_lags <- lapply(0:8, function(lag) {
  readRDS(paste0("output/san/seasonal/rnd_MSY_esc_annual_lag", lag, ".rds"))
})
names(stk_esc_lags) <- 0:8

stats_esc <- lapply(seq_along(stk_esc_lags), function(x) {
  #browser()
  qnts <- FLQuants("Catch" = seasonSums(catch(stk_esc_lags[[x]])[, ac(101:125)])/
                     refpts["msy", "yield"],
                   "Fbar" = seasonSums(fbar(stk_esc_lags[[x]])[, ac(101:125)])/
                     refpts["msy", "harvest"],
                   "SSB" = ssb(stk_esc_lags[[x]])[, ac(101:125),, 1]/
                     refpts["msy", "ssb"],
                   "TSB" = tsb(stk_esc_lags[[x]])[, ac(101:125),, 1]/
                     refpts["virgin", "biomass"],
                   "risk" = ssb(stk_esc_lags[[x]])[, ac(101:125),, 1] %=% NA_real_)
  qnts$risk[] <- as.numeric(ssb(stk_esc_lags[[x]])[, ac(101:125),, 1] < 163)
  qnts <- as.data.frame(qnts)
  qnts$lag <- as.numeric(names(stk_esc_lags)[x])
  return(qnts)
})
stats_esc <- do.call(rbind, stats_esc)

p_esc_catch <- stats_esc %>%
  filter(qname == "Catch") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(annual~catch/MSY), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 2.5)) +
  facet_grid(~ "Escapement strategy") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_esc_catch
p_esc_ssb <- stats_esc %>%
  filter(qname == "SSB") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(SSB/SSB[MSY]), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_esc_ssb
p_esc_tsb <- stats_esc %>%
  filter(qname == "TSB") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(fill = "grey", size = 0.2, show.legend = FALSE, adjust = 2) +
  #scale_fill_brewer(palette = "Set1") +
  geom_boxplot(fill = "white", width = 0.025, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(B/B[0]), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 2)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_esc_tsb
p_esc_risk <- stats_esc %>%
  filter(qname == "risk") %>%
  ggplot(aes(x = lag/4, y = data, group = lag)) +
  #geom_hline(yintercept = 1, colour = "grey") +
  #geom_violin(fill = "grey", size = 0.2, show.legend = FALSE) +
  geom_col(data = stats_esc %>%
             filter(qname == "risk") %>%
             group_by(lag) %>%
             summarise(data = mean(data)),
           aes(x = lag/4, y = data, group = lag), 
           fill = "grey", colour = "black", size = 0.2) +
  #scale_fill_brewer(palette = "Set1") +
  # geom_boxplot(fill = "white", width = 0.1, size = 0.2,
  #              outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
  #              outlier.fill = "transparent") +
  labs(y = expression(B[lim]~risk), x = "time lag [years]") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"))
p_esc_risk

p_esc_stats <- plot_grid(p_esc_catch, p_esc_ssb, p_esc_tsb, p_esc_risk, ncol = 1, 
                         align = "v", rel_heights = c(1, 1, 1, 1.2))
p_esc_stats
ggsave("output/plots/san_MSY_esc_lag_stats.png",
       width = 8.5, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_esc_lag_stats.pdf",
       width = 8.5, height = 10, units = "cm")

### combine HR and Escapement
# plot_grid(p_hr_stats, p_esc_stats, nrow = 1)
p_stats <- plot_grid(p_esc_catch, p_esc_ssb, p_esc_tsb, p_esc_risk,
                     p_catch, p_ssb, p_tsb, p_risk, 
                     ncol = 2, byrow = FALSE,
                     align = "v", rel_heights = c(1.17, 1, 1, 1.2))
ggsave("output/plots/san_MSY_hr_esc_lag_stats.png", plot = p_stats,
       width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_hr_esc_lag_stats.pdf", plot = p_stats,
       width = 17, height = 10, units = "cm")





