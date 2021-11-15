### ------------------------------------------------------------------------ ###
### find seasonal MSY harvest rate ####
### ------------------------------------------------------------------------ ###
source("funs.R")
library(FLasher)
library(doParallel)
library(tidyverse)

cl <- makeCluster(10)
registerDoParallel(cl)
. <- foreach(1:10) %dopar% {
  source("funs.R")
  library(FLasher)
}

### load deterministic OMs with 1 iteration
stk <- readRDS("input/san/stk_fmsy_1.rds")
sr <- readRDS("input/san/sr_1.rds")
residuals(sr)[,,, 1] <- 1 ### remove recruitment residuals

plot(stk)



# debugonce(mse_hr)
optimise_hr <- function(om_stk = stk, om_sr = sr, 
                        yrs = 101:200, seasons = 1:4, 
                        hr, lag = 0, catch_interval = 1,
                        verbose = FALSE, stat_yrs = 191:200,
                        return_all = FALSE) {
  res <- mse_hr(om_stk = om_stk, om_sr = om_sr, 
                yrs = yrs, seasons = seasons, 
                hr = hr, lag = lag, catch_interval = catch_interval,
                verbose = verbose)
  res_list <- list(TSB = median(tsb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
                   SSB = median(ssb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
                   Catch = median(apply(catch(res)[, ac(stat_yrs)], 2, sum), 
                                  na.rm = TRUE),
                   Fbar = median(apply(fbar(res)[, ac(stat_yrs)], 2, sum), 
                                 na.rm = TRUE),
                   Rec = median(rec(res)[, ac(stat_yrs),, 1], na.rm = TRUE))
  print(c(hr = hr, unlist(res_list)))
  if (isTRUE(return_all)) {
    res <- res_list
  } else {
    res <- mean(apply(catch(res)[, ac(stat_yrs)], 2, sum), 
                na.rm = TRUE)
    ### use mean so that if stock collapses in last year, objective value is
    ### reduced
  }
}

### check
if (FALSE) {
  out <- optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                     hr = 0.1, lag = 0, catch_interval = 1, verbose = TRUE, 
                     stat_yrs = 191:200, return_all = TRUE)
  out
}
### run a few values
hr_vals <- c(seq(0, 0.5, 0.02), 0.25, 0.27, 0.49)
hr_runs <- foreach(hr = hr_vals) %dopar% {
  optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 1, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE, hr = hr)
}
hr_runs <- as.data.frame(do.call(bind_rows, hr_runs))
hr_runs$hr <- hr_vals

hr_runs %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)

### MSY harvest rate
hr_res <- optimise(f = optimise_hr, 
                   om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                   lag = 0, catch_interval = 1, verbose = FALSE, 
                   stat_yrs = 191:200, return_all = FALSE,
                   interval = c(0, 0.5),
                   lower = 0, upper = 0.5,
                   maximum = TRUE,
                   tol = 0.00001)
saveRDS(hr_res, "input/san/hr_res_MSY.rds")
hr_res <- readRDS("input/san/hr_res_MSY.rds")
hr_runs %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = hr_res$maximum)
### run max catch
hr_res_MSY <- optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, 
                          seasons = 1:4, hr = hr_res$maximum, lag = 0, 
                          catch_interval = 1, verbose = TRUE,  
                          stat_yrs = 191:200, return_all = TRUE)
### add MSY run
hr_runs <- bind_rows(hr_runs,
                     c(unlist(hr_res_MSY), hr = hr_res$maximum))
saveRDS(hr_runs, "input/san/hr_runs.rds")
hr_runs <- readRDS("input/san/hr_runs.rds")
### run MSY projection and return stock
hr_res_MSY_stk <- mse_hr(om_stk = stk, om_sr = sr,
                         yrs = 101:200, seasons = 1:4,
                         hr = hr_res$maximum, lag = 0, catch_interval = 1,
                         verbose = TRUE)
plot(hr_res_MSY_stk)
saveRDS(hr_res_MSY_stk, "input/san/hr_res_MSY_stk.rds")
hr_res_MSY_stk <- readRDS("input/san/hr_res_MSY_stk.rds")


### do same with annual catch
### run a few values
hr_vals_annual <- seq(0, 0.3, 0.01)
hr_runs_annual <- foreach(hr = hr_vals_annual) %dopar% {
  optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 4, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE, hr = hr)
}
hr_runs_annual <- as.data.frame(do.call(bind_rows, hr_runs_annual))
hr_runs_annual$hr <- hr_vals_annual
hr_runs_annual %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)
### MSY harvest rate
hr_res_annual <- optimise(f = optimise_hr, 
                   om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                   lag = 0, catch_interval = 4, verbose = FALSE, 
                   stat_yrs = 191:200, return_all = FALSE,
                   interval = c(0, 0.22),
                   lower = 0, upper = 0.22,
                   maximum = TRUE,
                   tol = 0.00001)
saveRDS(hr_res_annual, "input/san/hr_res_MSY_annual.rds")
hr_res_annual <- readRDS("input/san/hr_res_MSY_annual.rds")
hr_runs_annual %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = hr_res_annual$maximum)
### run max catch
hr_res_annual_MSY <- optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, 
                          seasons = 1:4, hr = hr_res_annual$maximum, lag = 0, 
                          catch_interval = 4, verbose = TRUE,  
                          stat_yrs = 191:200, return_all = TRUE)
### add MSY run
hr_runs_annual <- bind_rows(hr_runs_annual,
                     c(unlist(hr_res_annual_MSY), hr = hr_res_annual$maximum))
saveRDS(hr_runs_annual, "input/san/hr_runs_annual.rds")
hr_runs_annual <- readRDS("input/san/hr_runs_annual.rds")
### run MSY projection and return stock
hr_res_MSY_annual_stk <- mse_hr(om_stk = stk, om_sr = sr,
                         yrs = 101:200, seasons = 1:4,
                         hr = hr_res_annual$maximum, lag = 0, catch_interval = 4,
                         verbose = TRUE)
plot(hr_res_MSY_annual_stk)
saveRDS(hr_res_MSY_annual_stk, "input/san/hr_res_MSY_annual_stk.rds")
hr_res_MSY_annual_stk <- readRDS("input/san/hr_res_MSY_annual_stk.rds")

### do same with biannual catch
### run a few values
hr_vals_biannual <- seq(0, 0.4, 0.02)
hr_runs_biannual <- foreach(hr = hr_vals_biannual) %dopar% {
  optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 2, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE, hr = hr)
}
hr_runs_biannual <- as.data.frame(do.call(bind_rows, hr_runs_biannual))
hr_runs_biannual$hr <- hr_vals_biannual
saveRDS(hr_runs_biannual, "input/san/hr_runs_biannual.rds")
hr_runs_biannual %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)
### MSY harvest rate
hr_res_biannual <- optimise(f = optimise_hr, 
                          om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                          lag = 0, catch_interval = 2, verbose = FALSE, 
                          stat_yrs = 191:200, return_all = FALSE,
                          interval = c(0, 0.4),
                          lower = 0, upper = 0.4,
                          maximum = TRUE,
                          tol = 0.00001)
saveRDS(hr_res_biannual, "input/san/hr_res_MSY_biannual.rds")
hr_res_biannual <- readRDS("input/san/hr_res_MSY_biannual.rds")
hr_runs_biannual %>% 
  ggplot(aes(x = hr, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = hr_res_biannual$maximum)
### run max catch
hr_res_biannual_MSY <- optimise_hr(om_stk = stk, om_sr = sr, yrs = 101:200, 
                                 seasons = 1:4, hr = hr_res_biannual$maximum, lag = 0, 
                                 catch_interval = 2, verbose = TRUE,  
                                 stat_yrs = 191:200, return_all = TRUE)
### add MSY run
hr_runs_biannual <- bind_rows(hr_runs_biannual,
                            c(unlist(hr_res_biannual_MSY), hr = hr_res_biannual$maximum))
saveRDS(hr_runs_biannual, "input/san/hr_runs_biannual.rds")
hr_runs_biannual <- readRDS("input/san/hr_runs_biannual.rds")
### run MSY projection and return stock
hr_res_MSY_biannual_stk <- mse_hr(om_stk = stk, om_sr = sr,
                                yrs = 101:200, seasons = 1:4,
                                hr = hr_res_biannual$maximum, lag = 0, catch_interval = 2,
                                verbose = TRUE)
plot(hr_res_MSY_biannual_stk)
saveRDS(hr_res_MSY_biannual_stk, "input/san/hr_res_MSY_biannual_stk.rds")
hr_res_MSY_biannual_stk <- readRDS("input/san/hr_res_MSY_biannual_stk.rds")



### combine quarterly/biannual/annual
bind_rows(hr_runs %>% mutate(step = "quarterly"),
          hr_runs_biannual %>% mutate(step = "biannual"),
          hr_runs_annual %>% mutate(step = "annual")) %>%
  mutate(step = factor(step, levels = c("quarterly", "biannual", "annual"))) %>%
  mutate(Rec = Rec/1000) %>%
  mutate(Fbar = ifelse(Fbar > 5, NA, Fbar)) %>%
  pivot_longer(c(SSB, Catch, Fbar, Rec)) %>%
  mutate(name = factor(name, levels = c("Rec", "SSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", 
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  ggplot(aes(x = hr, y = value, colour = step)) +
  geom_line(size = 0.4) +
  geom_vline(data = 
    data.frame(hr = c(hr_res$maximum, hr_res_biannual$maximum, 
                      hr_res_annual$maximum),
               step = c("quarterly", "biannual", "annual")),
    aes(xintercept = hr, colour = step), size = 0.3, linetype = "dotted",
    show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") + 
  # scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
  #                    name = "season") +
  # scale_linetype("target") +
  # coord_cartesian(ylim = c(0, NA), xlim = c(0, 0.75), expand = 1) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", nrow = 2) +
  labs(x = "Harvest rate [quarterly]") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        #legend.position = c(0.9, 0.3),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_HR_step_comparison.png",
       width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_HR_step_comparison.pdf",
       width = 17, height = 10, units = "cm")

### compare projections - quarterly, biannual, annual
qnts <- FLQuants(quarterly_Rec = rec(hr_res_MSY_stk),
                     quarterly_SSB = ssb(hr_res_MSY_stk),
                     quarterly_Catch = catch(hr_res_MSY_stk),
                     quarterly_Fbar = fbar(hr_res_MSY_stk),
                     biannual_Rec = rec(hr_res_MSY_biannual_stk),
                     biannual_SSB = ssb(hr_res_MSY_biannual_stk),
                     biannual_Catch = catch(hr_res_MSY_biannual_stk),
                     biannual_Fbar = fbar(hr_res_MSY_biannual_stk),
                     annual_Rec = rec(hr_res_MSY_annual_stk),
                     annual_SSB = ssb(hr_res_MSY_annual_stk),
                     annual_Catch = catch(hr_res_MSY_annual_stk),
                     annual_Fbar = fbar(hr_res_MSY_annual_stk))
as.data.frame(window(qnts, start = 149, end = 151)) %>%
  separate(col = qname, sep = "_", into = c("step", "quant")) %>%
  mutate(data = ifelse(quant == "Rec" & season != 1, 0, data)) %>%
  mutate(data = ifelse(quant == "Rec", data/1000, data),
         season = as.numeric(as.character(season)),
         year = as.numeric(as.character(year)),
         time = year + (season - 1)/4 - 150,
         step = factor(step, levels = c("quarterly", "biannual", "annual")),
         quant = factor(quant, levels = c("Rec", "SSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", 
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  # arrange(step, quant, time) %>%
  #mutate(time = year + (season - 1)/4 - 50) %>%
  ggplot(aes(x = time, y = data, colour = step, linetype = step)) +
  geom_line(size = 0.4) +
  scale_color_brewer(name = "step", palette = "Dark2") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
                     name = "season") +
  scale_linetype("step") +
  coord_cartesian(ylim = c(0, NA), xlim = c(0, 0.75), expand = 1) +
  facet_wrap(~ quant, scales = "free_y", strip.position = "left", nrow = 1) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.9, 0.3),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_hr_steps_projection.png",
       width = 17, height = 5, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_hr_steps_projection.pdf",
       width = 17, height = 5, units = "cm")



