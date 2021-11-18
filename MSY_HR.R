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

### check
if (FALSE) {
  out <- optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                     hr = 0.1, lag = 0, catch_interval = 1, verbose = TRUE, 
                     stat_yrs = 191:200, return_all = TRUE)
  out
}
### run a few values
hr_vals <- c(seq(0, 0.5, 0.02), 0.25, 0.27, 0.49)
hr_runs <- foreach(hr = hr_vals) %dopar% {
  optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
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
hr_res <- optimise(f = optimise_MP, 
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
hr_res_MSY <- optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, 
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
  optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
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
hr_res_annual <- optimise(f = optimise_MP, 
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
hr_res_annual_MSY <- optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, 
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
  optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
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
hr_res_biannual <- optimise(f = optimise_MP, 
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
hr_res_biannual_MSY <- optimise_MP(om_stk = stk, om_sr = sr, yrs = 101:200, 
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
  pivot_longer(c(SSB, TSB, Catch, Fbar, Rec)) %>%
  mutate(name = factor(name, levels = c("Rec", "SSB", "TSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  ggplot(aes(x = hr, y = value, colour = step)) +
  geom_line(size = 0.4) +
  geom_blank(data = 
               data.frame(name = factor(c("SSB [t]"), 
                                        levels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                                   "Catch [1000t]", "mean F (ages 1-3)")),
                          hr = 0, value = 1000, 
                          step = "annual")) +
  geom_vline(data = 
    data.frame(hr = c(hr_res$maximum, hr_res_biannual$maximum, 
                      hr_res_annual$maximum),
               step = c("quarterly", "biannual", "annual")),
    aes(xintercept = hr, colour = step), size = 0.3, linetype = "dashed",
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
        legend.position = c(0.8, 0.25),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_HR_step_comparison.png",
       width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_HR_step_comparison.pdf",
       width = 17, height = 10, units = "cm")

### compare projections - quarterly, biannual, annual
qnts <- FLQuants(quarterly_Rec = rec(hr_res_MSY_stk),
                 quarterly_SSB = ssb(hr_res_MSY_stk),
                 quarterly_TSB = tsb(hr_res_MSY_stk),
                 quarterly_Catch = catch(hr_res_MSY_stk),
                 quarterly_Fbar = fbar(hr_res_MSY_stk),
                 biannual_Rec = rec(hr_res_MSY_biannual_stk),
                 biannual_SSB = ssb(hr_res_MSY_biannual_stk),
                 biannual_TSB = tsb(hr_res_MSY_biannual_stk),
                 biannual_Catch = catch(hr_res_MSY_biannual_stk),
                 biannual_Fbar = fbar(hr_res_MSY_biannual_stk),
                 annual_Rec = rec(hr_res_MSY_annual_stk),
                 annual_SSB = ssb(hr_res_MSY_annual_stk),
                 annual_TSB = tsb(hr_res_MSY_annual_stk),
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
         quant = factor(quant, levels = c("Rec", "SSB", "TSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  # arrange(step, quant, time) %>%
  #mutate(time = year + (season - 1)/4 - 50) %>%
  ggplot(aes(x = time, y = data, colour = step)) +
  geom_line(size = 0.4) +
  scale_color_brewer(name = "step", palette = "Dark2") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
                     name = "season") +
  scale_linetype("step") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-0.085, 0.835), expand = 1) +
  facet_wrap(~ quant, scales = "free_y", strip.position = "left", nrow = 1) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.08, 0.8),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_hr_steps_projection.png",
       width = 17, height = 5, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_hr_steps_projection.pdf",
       width = 17, height = 5, units = "cm")


### ------------------------------------------------------------------------ ###
### escapement strategy: find max long-term catch ####
### ------------------------------------------------------------------------ ###

# debugonce(mse_hr)


### run a few values
esc_vals <- c(seq(0, 900, 50), seq(1000, 2000, 100))
esc_runs <- foreach(esc_biomass = esc_vals) %dopar% {
  optimise_MP(MP = "escapement", 
              om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 1, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE,
              target = esc_biomass)
}
esc_runs <- as.data.frame(do.call(bind_rows, esc_runs))
esc_runs$esc_biomass <- esc_vals

esc_runs %>% 
  ggplot(aes(x = esc_biomass, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)

### MSY escapement
trace_env <- new.env()
#source("funs.R")
#get("res_trace", envir = trace_env)
#debugonce(optimise_MP)
esc_res <- optimise(f = optimise_MP, MP = "escapement", 
                    om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                    lag = 0, catch_interval = 1, verbose = FALSE, 
                    stat_yrs = 191:200, return_all = FALSE,
                    trace = TRUE, trace_env = trace_env,
                    interval = c(0, 1000),
                    lower = 0, upper = 1000,
                    maximum = TRUE,
                    tol = 0.1)
saveRDS(esc_res, "input/san/esc_res_MSY.rds")
esc_res <- readRDS("input/san/esc_res_MSY.rds")
### add MSY run
esc_runs <- bind_rows(esc_runs,
                      bind_rows(get("res_trace", envir = trace_env)))
saveRDS(esc_runs, "input/san/esc_runs.rds")
esc_runs <- readRDS("input/san/esc_runs.rds")
esc_runs %>% 
  ggplot(aes(x = target, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = esc_res$maximum)
### run MSY projection and return stock
esc_res_MSY_stk <- mse_loop(MP = "escapement",
                         om_stk = stk, om_sr = sr,
                         yrs = 101:200, seasons = 1:4,
                         target = esc_res$maximum, 
                         lag = 0, catch_interval = 1,
                         verbose = TRUE)
plot(esc_res_MSY_stk)
saveRDS(esc_res_MSY_stk, "input/san/esc_res_MSY_stk.rds")
esc_res_MSY_stk <- readRDS("input/san/esc_res_MSY_stk.rds")


### same with biannual catch
### run a few values
esc_vals_biannual <- c(seq(0, 950, 50), seq(1000, 2000, 100))
esc_runs_biannual <- foreach(esc_biomass = esc_vals_biannual) %dopar% {
  optimise_MP(MP = "escapement", 
              om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 2, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE,
              target = esc_biomass)
}
esc_runs_biannual <- as.data.frame(do.call(bind_rows, esc_runs_biannual))
esc_runs_biannual %>% 
  ggplot(aes(x = target, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)

### MSY escapement
trace_env <- new.env()
#source("funs.R")
#get("res_trace", envir = trace_env)
#debugonce(optimise_MP)
esc_res_biannual <- optimise(f = optimise_MP, MP = "escapement", 
                    om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                    lag = 0, catch_interval = 2, verbose = FALSE, 
                    stat_yrs = 191:200, return_all = FALSE,
                    trace = TRUE, trace_env = trace_env, 
                    objective_stat = median,
                    interval = c(0, 1000),
                    lower = 0, upper = 1000,
                    maximum = TRUE,
                    tol = 0.1)
saveRDS(esc_res_biannual, "input/san/esc_res_biannual_MSY.rds")
esc_res_biannual <- readRDS("input/san/esc_res_biannual_MSY.rds")
### add MSY run
esc_runs_biannual <- bind_rows(esc_runs_biannual,
                      bind_rows(get("res_trace", envir = trace_env)))
saveRDS(esc_runs_biannual, "input/san/esc_runs_biannual.rds")
esc_runs_biannual <- readRDS("input/san/esc_runs_biannual.rds")
esc_runs_biannual %>% 
  ggplot(aes(x = target, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = esc_res_biannual$maximum)
### run MSY projection and return stock
esc_res_MSY_stk_biannual <- mse_loop(MP = "escapement",
                                     om_stk = stk, om_sr = sr,
                                     yrs = 101:200, seasons = 1:4,
                                     target = esc_res_biannual$maximum, 
                                     lag = 0, catch_interval = 2,
                                     verbose = TRUE)
plot(esc_res_MSY_stk_biannual)
saveRDS(esc_res_MSY_stk_biannual, "input/san/esc_res_MSY_stk_biannual.rds")
esc_res_MSY_stk_biannual <- readRDS("input/san/esc_res_MSY_stk_biannual.rds")



### same with annual catch
### run a few values
esc_vals_annual <- unique(c(seq(0, 1000, 50), seq(1000, 2000, 100)))
esc_runs_annual <- foreach(esc_biomass = esc_vals_annual) %dopar% {
  optimise_MP(MP = "escapement", 
              om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
              lag = 0, catch_interval = 4, verbose = FALSE, 
              stat_yrs = 191:200, return_all = TRUE,
              target = esc_biomass, objective_stat = median)
}
esc_runs_annual <- as.data.frame(do.call(bind_rows, esc_runs_annual))
esc_runs_annual %>% 
  ggplot(aes(x = target, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8)
### MSY escapement
trace_env <- new.env()
esc_res_annual <- optimise(f = optimise_MP, MP = "escapement", 
                             om_stk = stk, om_sr = sr, yrs = 101:200, seasons = 1:4, 
                             lag = 0, catch_interval = 4, verbose = FALSE, 
                             stat_yrs = 191:200, return_all = FALSE,
                             trace = TRUE, trace_env = trace_env, 
                             objective_stat = median,
                             interval = c(0, 1000),
                             lower = 0, upper = 1000,
                             maximum = TRUE,
                             tol = 0.1)
saveRDS(esc_res_annual, "input/san/esc_res_annual_MSY.rds")
esc_res_annual <- readRDS("input/san/esc_res_annual_MSY.rds")
### add MSY run
esc_runs_annual <- bind_rows(esc_runs_annual,
                               bind_rows(get("res_trace", envir = trace_env)))
saveRDS(esc_runs_annual, "input/san/esc_runs_annual.rds")
esc_runs_annual <- readRDS("input/san/esc_runs_annual.rds")
esc_runs_annual %>% 
  ggplot(aes(x = target, y = Catch)) +
  geom_line() +
  theme_bw(base_size = 8) +
  geom_vline(xintercept = esc_res_annual$maximum)
### run MSY projection and return stock
esc_res_MSY_stk_annual <- mse_loop(MP = "escapement",
                                     om_stk = stk, om_sr = sr,
                                     yrs = 101:200, seasons = 1:4,
                                     target = esc_res_annual$maximum, 
                                     lag = 0, catch_interval = 4,
                                     verbose = TRUE)
plot(esc_res_MSY_stk_annual)
saveRDS(esc_res_MSY_stk_annual, "input/san/esc_res_MSY_stk_annual.rds")
esc_res_MSY_stk_annual <- readRDS("input/san/esc_res_MSY_stk_annual.rds")


### combine quarterly/biannual/annual
bind_rows(esc_runs %>% mutate(step = "quarterly"),
          esc_runs_biannual %>% mutate(step = "biannual"),
          esc_runs_annual %>% mutate(step = "annual")) %>%
  mutate(step = factor(step, levels = c("quarterly", "biannual", "annual"))) %>%
  mutate(Rec = Rec/1000) %>%
  mutate(Fbar = ifelse(Fbar > 5, NA, Fbar)) %>%
  pivot_longer(c(TSB, SSB, Catch, Fbar, Rec)) %>%
  mutate(name = factor(name, levels = c("Rec", "SSB", "TSB", "Catch", "Fbar"),
                       labels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                  "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  ggplot(aes(x = target, y = value, colour = step)) +
  geom_line(size = 0.4) +
  geom_blank(data = 
    data.frame(name = factor(c("SSB [t]"), 
                             levels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                        "Catch [1000t]", "mean F (ages 1-3)")),
                               target = 0, value = 1000, 
                               step = "annual")) +
  geom_vline(data = 
               data.frame(target = c(esc_res$maximum, esc_res_biannual$maximum, 
                                 esc_res_annual$maximum),
                          step = c("quarterly", "biannual", "annual")),
             aes(xintercept = target, colour = step), size = 0.3, 
             linetype = "dashed",
             show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") + 
  # scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
  #                    name = "season") +
  # scale_linetype("target") +
  # coord_cartesian(ylim = c(0, NA), xlim = c(0, 0.75), expand = 1) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left", nrow = 2) +
  labs(x = "Escapement biomass [t]") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.8, 0.25),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_esc_step_comparison.png",
       width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_esc_step_comparison.pdf",
       width = 17, height = 10, units = "cm")

### compare projections - quarterly, biannual, annual
qnts <- FLQuants(quarterly_Rec = rec(esc_res_MSY_stk),
                 quarterly_SSB = ssb(esc_res_MSY_stk),
                 quarterly_TSB = tsb(esc_res_MSY_stk),
                 quarterly_Catch = catch(esc_res_MSY_stk),
                 quarterly_Fbar = fbar(esc_res_MSY_stk),
                 biannual_Rec = rec(esc_res_MSY_stk_biannual),
                 biannual_SSB = ssb(esc_res_MSY_stk_biannual),
                 biannual_TSB = tsb(esc_res_MSY_stk_biannual),
                 biannual_Catch = catch(esc_res_MSY_stk_biannual),
                 biannual_Fbar = fbar(esc_res_MSY_stk_biannual),
                 annual_Rec = rec(esc_res_MSY_stk_annual),
                 annual_SSB = ssb(esc_res_MSY_stk_annual),
                 annual_TSB = tsb(esc_res_MSY_stk_annual),
                 annual_Catch = catch(esc_res_MSY_stk_annual),
                 annual_Fbar = fbar(esc_res_MSY_stk_annual))
as.data.frame(window(qnts, start = 149, end = 151)) %>%
  separate(col = qname, sep = "_", into = c("step", "quant")) %>%
  mutate(data = ifelse(quant == "Rec" & season != 1, 0, data)) %>%
  mutate(data = ifelse(quant == "Rec", data/1000, data),
         season = as.numeric(as.character(season)),
         year = as.numeric(as.character(year)),
         time = year + (season - 1)/4 - 150,
         step = factor(step, levels = c("quarterly", "biannual", "annual")),
         quant = factor(quant, levels = c("Rec", "SSB", "TSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", "TSB [t]",
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  # arrange(step, quant, time) %>%
  #mutate(time = year + (season - 1)/4 - 50) %>%
  ggplot(aes(x = time, y = data, colour = step)) +
  geom_line(size = 0.4) +
  scale_color_brewer(name = "step", palette = "Dark2") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
                     name = "season") +
  scale_linetype("step") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-0.085, 0.835), expand = 1) +
  facet_wrap(~ quant, scales = "free_y", strip.position = "left", nrow = 1) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.08, 0.8),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave("output/plots/san_MSY_esc_steps_projection.png",
       width = 17, height = 5, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_esc_steps_projection.pdf",
       width = 17, height = 5, units = "cm")

