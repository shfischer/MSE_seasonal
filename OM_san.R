### ------------------------------------------------------------------------ ###
### section ####
### ------------------------------------------------------------------------ ###

library(FLCore)
library(FLBRP)
library(FLife)
library(ggplot2)
library(dplyr)
library(tidyr)
library(FLasher)
library(foreach)
# remotes::install_github("shfischer/FLasher", ref = "seasons_tmp",
#                         INSTALL_opts = "--no-multiarch")
### use shfischer/FLasher branch seasons_tmp
source("funs.R")

### ------------------------------------------------------------------------ ###
### load life-history parameters ####
### ------------------------------------------------------------------------ ###
stocks <- read.csv(file = "input/stocks.csv")
stocks <- read.csv(file = "../../data-limited/MSE_seasonal/input/stocks.csv")
stocks$t0[is.na(stocks$t0)] <- -0.1
stocks$s <- 0.75

### ------------------------------------------------------------------------ ###
### sandeel: 4 seasons ####
### ------------------------------------------------------------------------ ###

### example stock: sandeel
stocks[stocks$stock == "san", ]
lh <- FLPar(a = stocks[stocks$stock == "san", "a"], 
            b = stocks[stocks$stock == "san", "b"], 
            linf = stocks[stocks$stock == "san", "linf"], 
            l50 = stocks[stocks$stock == "san", "l50"], 
            k = stocks[stocks$stock == "san", "k"], 
            t0 = stocks[stocks$stock == "san", "t0"], 
            s = stocks[stocks$stock == "san", "s"])
lh <- lhPar(lh)
max_age <- ceiling(log(0.05)/(-c(lh["k"])) + c(lh["t0"]))

# brp <- lhEql(params = lh)
# plot(brp)

ages <- seq(0, 5, 0.0001)
lengths <- data.frame(age = ages,
                      length = vonB(age = ages, params = lh))
lengths <- bind_rows(lengths %>% 
                       filter(age %in% seq(0, 5, 1)) %>%
                       mutate(step = "annual"),
                     lengths %>% 
                       filter(age %in% seq(0, 5, 0.25)) %>%
                       mutate(step = "seasonal"),
                     lengths %>% 
                       mutate(step = "continuous")) %>%
  mutate(step = factor(step, levels = c("continuous", "seasonal", "annual")))
#lengths %>%
ggplot() +
  #geom_line(data = lengths, aes(x = age, y = length)) +
  geom_step(data = lengths,
            aes(x = age, y = length, colour = step), size = 0.4) +
  scale_colour_brewer(palette = "Dark2") +
  #scale_linetype_manual(values = c("dotted", "solid", "2121")) +
  ylim(c(0, NA)) +
  labs(x = "age [years]", y = "length [cm]") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.8, 0.5),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank())
# ggsave("output/plots/san_length.png",
#        width = 8.5, height = 6, units = "cm", dpi = 300, type = "cairo")
# ggsave("output/plots/san_length.pdf",
#        width = 8.5, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### create seasonal OM (deterministic) ####
### ------------------------------------------------------------------------ ###
nseasons <- 4
range <- list(min = 0, max = 5, 
              minfbar = 1, maxfbar = 3,
              plusgrup = 5,
              seasons = nseasons)
params <- lh

# debugonce(lhEql, signature = "FLPar")
# brp <- lhEql(params = lh, range = range)

### timing
spwn <- 0 ### spawning at beginning of year
### fishing occurs throughout
### -> use mid-year (or mid-season) lengths/weights 
fish <- 0.5 
midyear <- 0.5

### ages as integer vector
ages_int <- seq(from = floor(range$min), to = ceiling(range$max))
### FLQuant template with all dimensions
flq <- FLQuant(NA, dimnames = list(age = ages_int,
                                   season = seq(range$seasons)))
### all ages (including seasons) as FLQuant
ages_full <- seq(from = min(ages_int), 
                to = max(ages_int) + 1 - 1/range$seasons, 
                by = 1/range$seasons)
ages_full <- flq %=% c(matrix(ages_full, ncol = range$seasons, byrow = TRUE))

### spawning period - proportion of F/M per season before spawning
season_times <- seq(range$seasons)/range$seasons - 1/range$seasons
m.spwn <- sapply(season_times, function(x) {
  ifelse(spwn %in% season_times,
         ifelse(x < spwn, 1, 0),
         ifelse(x > spwn, 0,
                ifelse((spwn - x) / (1/range$seasons) > 1, 1,
                ((spwn %% (1/range$seasons))/(1/range$seasons)))))
})
m.spwn <- flq %=% rep(m.spwn, each = dim(flq)[1])
harvest.spwn <- m.spwn

### lengths
lengths_stk <- FLife::vonB(ages_full, params = params)
lengths_catch <- FLife::vonB(ages_full + fish/nseasons, params = params)
lengths_mid <- FLife::vonB(ages_full + 1/range$seasons*midyear, params = params)

### weights
weights_stk <- FLife::len2wt(length = lengths_stk, params = params)
weights_catch <- FLife::len2wt(length = lengths_catch, params = params)
weights_mid <- FLife::len2wt(length = lengths_mid, params = params)

### maturity
### use value at beginning of year/season because this is when SSB is calculated
mat. <- FLife::logistic(age = ages_full, params = params)
if (dims(mat.)["min"] == 0) mat.[1,,, 1] <- 0 ### no self-spawning

### fisheries selectivity
sel. <- FLife::dnormal(age = ages_full + fish/nseasons, params = params)
if (dims(sel.)["min"] == 0) sel.[1,,, 1] <- 0 ### recruits not fished
units(sel.) <- "f"

### natural mortality
### define Gislason natural mortality function
### FLife::gislason is wrong
gisM <- function(length, params) {
  exp(params["m1"] %+% (params["m2"] %*% log(length)) %+% 
        (params["m3"] %*% log(params["linf"])) %+% log(params["k"]))
}
### /seasons because M is annual value
m. <- gisM(length = lengths_mid, params = params)/range$seasons

### create FLBRP template
brp <- FLBRP(stock.wt        = weights_stk,
             landings.wt     = weights_catch,
             discards.wt     = weights_catch,
             bycatch.wt      = weights_catch,
             mat             = mat.,
             m               = m.,
             landings.sel    = sel.,
             discards.sel    = sel. %=% 0,
             bycatch.harvest = sel. %=% 0,
             harvest.spwn    = harvest.spwn,
             m.spwn          = m.spwn,
             availability    = weights_stk %=% 1,
             range           = unlist(range))

### set up unfished stock
# debugonce(N3)
### all metrics based on SSB=1000 when F=0
SSB0 <- 1000
Ri <- 1e+06
Ni <- N(object = flq, m = m., rec = Ri, eq_years = 100)
### SSB in 1st season (spawning time)
SSBi <- sum((Ni * mat. * weights_stk)[,,, 1]) 
spr0 <- SSBi/Ri
R0 <- SSB0 * 1/spr0

### unfished stock
N0 <- N(object = flq, m = m., rec = R0, eq_years = 100)
sum((N0 * mat. * weights_stk)[,,, 1])

### define recruitment model
model(brp) <- bevholt()$model
a <- (4 * params["s"] * R0)/(5 * params["s"] - 1)
b <- (SSB0 * (1 - params["s"]))/(5 * params["s"] - 1)
sr <- as(brp, "FLSR")
params(sr) <- FLPar(NA, dimnames = list(params = c("a","b"),
                                        season = seq(4), iter = 1))
params(sr)[,1] <- c(a, b)

### prepare OM stock
stk <- as(brp, "FLStock")[, 1:101]
for (y in 1) stock.n(stk)[, y] <- N0 # do not fill later years -> weird things happen
stock(stk) <- computeStock(stk)
landings.n(stk) <- discards.n(stk) <- catch.n(stk) <- 0
landings.n(stk)[, 2:101] <- 1e-6
catch.n(stk)[, 2:101] <- 1e-6
catch(stk) <- computeCatch(stk)
landings(stk) <- computeLandings(stk)
discards(stk) <- computeDiscards(stk)
harvest(stk)[, 1] <- 0
for (y in 2:101) harvest(stk)[, y] <- sel.

saveRDS(stk, "input/san/stk_seasonal_deterministic.rds")
saveRDS(sr, "input/san/sr_seasonal_deterministic.rds")

### ------------------------------------------------------------------------ ###
### projections ####
### ------------------------------------------------------------------------ ###

# debugonce(fwd, signature = c("FLStock", "missing", "fwdControl"))

### zero fishing
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 0.0))
stk_fwd0 <- fwd(stk, control = ctrl, sr = sr, effort_max = 1e+9)
### fish 0.4
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 0.4/4))
stk_fwd0.1 <- fwd(stk, control = ctrl, sr = sr)
### fish 1
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 1/4))
stk_fwd1 <- fwd(stk, control = ctrl, sr = sr)
### fish 4
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 4/4))
stk_fwd5 <- fwd(stk, control = ctrl, sr = sr, effort_max = 1e+9)
### fish 10
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 10))
stk_fwd_ <- fwd(stk, control = ctrl, sr = sr, effort_max = 1e+9)

# plot(FLStocks("F=0" = stk_fwd0, "F=0.4" = stk_fwd0.1, "F=1" = stk_fwd1,
#               "F=4" = stk_fwd5)) +
#   theme_bw()
plot(window(FLStocks("F=0" = stk_fwd0, "F=0.4" = stk_fwd0.1, "F=1" = stk_fwd1,
                     "F=4" = stk_fwd5), end = 10)) +
  theme_bw()
# ggsave("output/plots/san_check_seasonal_F.png",
#        width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
# ggsave("output/plots/san_check_seasonal_F.pdf",
#        width = 17, height = 10, units = "cm")
plot(window(FLStocks("F=0" = stk_fwd0, "F=0.4" = stk_fwd0.1, "F=1" = stk_fwd1,
                     "F=4" = stk_fwd5), start = 95, end = 101)) +
  theme_bw()

### ------------------------------------------------------------------------ ###
### find MSY ####
### ------------------------------------------------------------------------ ###

fishF <- function(stk, sr, quant = c("fbar"), value, 
                  return_all = FALSE, maxF = 5) {
  # browser()
  ctrl_tmp <- fwdControl(data.frame(year = rep(2:101, each = 4),
                                    season = 1:4,
                                    quant = quant,
                                    value = value/dim(stk)[4]))
  stk_fwd_tmp <- fwd(stk, control = ctrl_tmp, sr = sr, effort_max = 1e+9,
                     maxF = maxF)
  # plot(stk_fwd_tmp)
  if (isTRUE(return_all)) {
    res <- list(TSB = median(tsb(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE),
                SSB = median(ssb(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE),
                Catch = median(apply(catch(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                               na.rm = TRUE),
                Fbar = median(apply(fbar(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                              na.rm = TRUE),
                Rec = median(rec(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE))
  } else {
    res <- mean(apply(catch(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                  na.rm = TRUE)
    ### use mean so that if stock collapses in last year, objective value is
    ### reduced
  }
}
fishF(stk = stk, sr = sr, quant = "fbar", value = 0, return_all = TRUE)
fishF(stk = stk, sr = sr, quant = "catch", value = 0, return_all = TRUE)
### run some values
runs <- lapply(X = seq(0, 5, 0.1), FUN = fishF, stk = stk, sr = sr, return_all = TRUE)
runs <- do.call(rbind, runs)
runs <- as.data.frame(runs)
# plot(runs$Fbar, runs$SSB, type = "l")
# plot(runs$Fbar, runs$Catch, type = "l")
# plot(runs$Fbar, runs$Rec, type = "l")
# plot(runs$SSB, runs$Rec, type = "l")
for (i in 1:ncol(runs)) runs[, i] <- unlist(runs[, i])
saveRDS(runs, file = "output/san_MSY_runs.rds")
runs <- readRDS("output/san_MSY_runs.rds")

### get MSY with seasonal F target
res <- optimise(f = fishF, stk = stk, sr = sr, quant = "fbar", return_all = FALSE,
                interval = c(0, 5),
                lower = 0, upper = 4,
                maximum = TRUE,
                tol = 0.00001)
saveRDS(res, file = "output/san_MSY_res.rds")
res <- readRDS("output/san_MSY_res.rds")
MSY_vals <- fishF(stk = stk, sr = sr, val = res$maximum, return_all = TRUE)
# $TSB
# [1] 983.4669
# $SSB
# [1] 221.2179
# $Catch
# [1] 457.1125
# $Fbar
# [1] 1.600873
# $Rec
# [1] 14414.4

### projection
ctrl_MSY <- fwdControl(data.frame(year = rep(2:101, each = 4),
                                       season = 1:4,
                                       quant = "fbar",
                                       value = res$maximum/4))
stk_fwd_MSY <- fwd(stk, control = ctrl_MSY, sr = sr, 
                   effort_max = 1e+9, maxF = 5)
plot(stk_fwd_MSY)

refpts <- FLPar(NA, dimnames = list(params = c("virgin","msy", "crash"),
                                    refpt = c("harvest", "yield", "rec", "ssb", 
                                              "biomass"),
                                    iter = 1))
refpts["virgin", ] <- unlist(runs[1, c("Fbar", "Catch", "Rec", "SSB", "TSB")])
refpts["msy", ] <- unlist(MSY_vals[c("Fbar", "Catch", "Rec", "SSB", "TSB")])
refpts["crash", ] <- unlist(runs[head(which(runs$SSB < 1), 1), 
                                 c("Fbar", "Catch", "Rec", "SSB", "TSB")])
saveRDS(refpts, file = "input/san/refpts.rds")
refpts <- readRDS("input/san/refpts.rds")

runs %>%
  pivot_longer(c(SSB, Catch, Rec)) %>%
  mutate(name = factor(name, levels = c("Catch", "SSB", "Rec"))) %>%
  ggplot(aes(x = Fbar, y = value)) +
  geom_line(size = 0.4) +
  geom_vline(xintercept = c(refpts["msy", "harvest"]), 
             linetype = "dashed", size = 0.3) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left") +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = TRUE) +
  geom_blank(data = data.frame(Fbar = 1,
                               value = c(420, 1000, 16500),
                               name = factor(c("Catch", "SSB", "Rec")))) +
  labs(x = "mean F (ages 1-3)") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 8))
ggsave("output/plots/san_seasonal_MSY.png",
       width = 17, height = 6, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_seasonal_MSY.pdf",
       width = 17, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### find MSY - with same catch in all seasons ####
### ------------------------------------------------------------------------ ###

### get MSY with seasonal F target
res_const <- optimise(f = fishF, stk = stk, sr = sr, quant = "catch", 
                      return_all = FALSE,
                      interval = c(0, 1000),
                      lower = 0, upper = 1000,
                      maximum = TRUE,
                      tol = 0.00001)
saveRDS(res_const, file = "output/san_MSY_cont_catch_res.rds")
res_const <- readRDS("output/san_MSY_cont_catch_res.rds")

### stats
MSY_vals_const <- fishF(stk = stk, sr = sr, quant = "catch", 
                         value = res_const$maximum, return_all = TRUE)
### projection
ctrl_MSY_const <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "catch",
                              value = 460/4))
stk_fwd_MSY_const <- fwd(stk, control = ctrl_MSY_const, sr = sr, 
                         effort_max = 1e+9, maxF = 5)
plot(stk_fwd_MSY_const)
plot(window(stk_fwd_MSY_const, start = 50, end = 55))
plot(window(tsb(stk_fwd_MSY_const), start = 50, end = 55))

### compare Fmsy (constant through year) and MSY (constant through year)
qnts_MSY <- FLQuants(F_Rec = rec(stk_fwd_MSY),
                     F_SSB = ssb(stk_fwd_MSY),
                     F_Catch = catch(stk_fwd_MSY),
                     F_Fbar = fbar(stk_fwd_MSY),
                     C_Rec = rec(stk_fwd_MSY_const),
                     C_SSB = ssb(stk_fwd_MSY_const),
                     C_Catch = catch(stk_fwd_MSY_const),
                     C_Fbar = fbar(stk_fwd_MSY_const))
as.data.frame(window(qnts_MSY, start = 49, end = 51)) %>%
  separate(col = qname, sep = "_", into = c("target", "quant")) %>%
  mutate(data = ifelse(quant == "Rec" & season != 1, 0, data)) %>%
  mutate(data = ifelse(quant == "Rec", data/1000, data),
         season = as.numeric(as.character(season)),
         year = as.numeric(as.character(year)),
         time = year + (season - 1)/4 - 50,
         target = factor(target, levels = c("F", "C"), 
                         labels = c("F", "Catch")),
         quant = factor(quant, levels = c("Rec", "SSB", "Catch", "Fbar"),
                        labels = c("Recruits [1000s]", "SSB [t]", 
                                   "Catch [1000t]", "mean F (ages 1-3)"))) %>%
  # arrange(target, quant, time) %>%
  #mutate(time = year + (season - 1)/4 - 50) %>%
  ggplot(aes(x = time, y = data, colour = target, linetype = target)) +
  geom_line(size = 0.4) +
  scale_color_brewer(name = "target", palette = "Dark2") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c(1, 2, 3, 4),
                     name = "season") +
  scale_linetype("target") +
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
ggsave("output/plots/san_MSY_seasonal_F_vs_catch.png",
       width = 17, height = 5, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY_seasonal_F_vs_catch.pdf",
       width = 17, height = 5, units = "cm")


### ------------------------------------------------------------------------ ###
### expand stock ####
### ------------------------------------------------------------------------ ###
# saveRDS(stk, "input/san/stk.rds")
stk <- readRDS("input/san/stk.rds")
# saveRDS(sr, "input/san/sr.rds")
sr <- readRDS("input/san/sr.rds")
n_iter <- 500

stk_fwd <- propagate(stf(stk, 99), n_iter)
sr_fwd <- sr

### ------------------------------------------------------------------------ ###
### recruitment residuals ####
### ------------------------------------------------------------------------ ###

### create recruitment residuals for (historical) projection
set.seed(0)
residuals(sr_fwd) <- rec(stk_fwd) %=% NA_real_
residuals(sr_fwd)[,,, 1] <- rlnoise(dim(stk_fwd)[6], rec(stk_fwd)[,,, 1] %=% 0, 
                                    sd = 0.6, b = 0)
saveRDS(sr_fwd, file = paste0("input/san/sr_", n_iter, ".rds"))

### ------------------------------------------------------------------------ ###
### create fishing histories ####
### ------------------------------------------------------------------------ ###
refpts <- readRDS("input/san/refpts.rds")

### fish at F=Fmsy
ctrl_fmsy <- FLQuant(c(refpts["msy", "harvest"])/4,
                     dimnames = list(year = 2:100, season = 1:4))
ctrl_fmsy <- as(FLQuants(fbar = ctrl_fmsy), "fwdControl")
stk_fwd_fmsy <- fwd(stk_fwd, control = ctrl_fmsy, sr = sr_fwd, 
                  deviances = residuals(sr_fwd),# %=% 1,
                  effort_max = 1e+9)
plot(window(stk_fwd_fmsy, end = 100))
plot(window(stk_fwd_fmsy, start = 90, end = 100))
saveRDS(stk_fwd_fmsy, file = paste0("input/san/stk_fmsy_", n_iter, ".rds"))
stk_fwd_fmsy <- readRDS(paste0("input/san/stk_fmsy_", n_iter, ".rds"))


### one-way fishing history
ctrl_ow <- fhist_one_way(n_iter = n_iter, yrs_hist = 100,
                         yrs_const = 75,
                         yrs_increase = 25, 
                         f0 = c(refpts["msy", "harvest"]) * 0.5,
                         fmax = c(refpts["crash", "harvest"]) * 0.8)
ctrl_ow <- ctrl_ow/4 ### divide by seasons because F is annual
ctrl_ow <- FLCore::expand(ctrl_ow, season = 1:4)
ctrl_ow <- as(FLQuants(fbar = ctrl_ow[, -1]), "fwdControl")
stk_fwd_ow <- fwd(stk_fwd, control = ctrl_ow, sr = sr_fwd, 
                  deviances = residuals(sr_fwd),# %=% 1,
                  effort_max = 1e+9)
plot(window(stk_fwd_ow, end = 100))
plot(window(stk_fwd_ow, start = 90, end = 100))
saveRDS(stk_fwd_ow, file = paste0("input/san/stk_ow_", n_iter, ".rds"))
stk_fwd_ow <- readRDS(paste0("input/san/stk_ow_", n_iter, ".rds"))

### random fishing history
set.seed(2)
ctrl_rnd <- fhist_random(n_iter = n_iter, yrs_hist = 100, 
                         min = 0,
                         max = c(refpts["crash", "harvest"]))
ctrl_rnd <- ctrl_rnd/4 ### divide by seasons because F is annual
ctrl_rnd <- FLCore::expand(ctrl_rnd, season = 1:4)
ctrl_rnd <- as(FLQuants(fbar = ctrl_rnd[, -1]), "fwdControl")
stk_fwd_rnd <- fwd(stk_fwd, control = ctrl_rnd, sr = sr_fwd, 
                   deviances = residuals(sr_fwd),
                   effort_max = 1e+9)
plot(window(stk_fwd_rnd, end = 100), iter = 1:5)
plot(window(stk_fwd_rnd, start = 90, end = 100), iter = 1:5)
saveRDS(stk_fwd_rnd, file = paste0("input/san/stk_rnd_", n_iter, ".rds"))
stk_fwd_rnd <- readRDS(paste0("input/san/stk_rnd_", n_iter, ".rds"))



### ------------------------------------------------------------------------ ###
### create OM again but without seasonal growth ####
### ------------------------------------------------------------------------ ###
nseasons <- 1
range <- list(min = 0, max = 5, 
              minfbar = 1, maxfbar = 3,
              plusgrup = 5,
              seasons = nseasons)
params <- lh

# debugonce(lhEql, signature = "FLPar")
# brp <- lhEql(params = lh, range = range)

### timing
spwn <- 0 ### spawning at beginning of year
### fishing occurs throughout
### -> use mid-year (or mid-season) lengths/weights 
fish <- 0.5 
midyear <- 0.5

### ages as integer vector
ages_int <- seq(from = floor(range$min), to = ceiling(range$max))
### FLQuant template with all dimensions
flq <- FLQuant(NA, dimnames = list(age = ages_int,
                                   season = seq(range$seasons)))
### all ages (including seasons) as FLQuant
ages_full <- seq(from = min(ages_int), 
                 to = max(ages_int) + 1 - 1/range$seasons, 
                 by = 1/range$seasons)
ages_full <- flq %=% c(matrix(ages_full, ncol = range$seasons, byrow = TRUE))

### spawning period - proportion of F/M per season before spawning
season_times <- seq(range$seasons)/range$seasons - 1/range$seasons
m.spwn <- sapply(season_times, function(x) {
  ifelse(spwn %in% season_times,
         ifelse(x < spwn, 1, 0),
         ifelse(x > spwn, 0,
                ifelse((spwn - x) / (1/range$seasons) > 1, 1,
                       ((spwn %% (1/range$seasons))/(1/range$seasons)))))
})
m.spwn <- flq %=% rep(m.spwn, each = dim(flq)[1])
harvest.spwn <- m.spwn

### lengths
lengths_stk <- FLife::vonB(ages_full, params = params)
lengths_catch <- FLife::vonB(ages_full + fish/nseasons, params = params)
lengths_mid <- FLife::vonB(ages_full + 1/range$seasons*midyear, params = params)

### weights
weights_stk <- FLife::len2wt(length = lengths_stk, params = params)
weights_catch <- FLife::len2wt(length = lengths_catch, params = params)
weights_mid <- FLife::len2wt(length = lengths_mid, params = params)

### maturity
### use value at beginning of year/season because this is when SSB is calculated
mat. <- FLife::logistic(age = ages_full, params = params)
if (dims(mat.)["min"] == 0) mat.[1,,, 1] <- 0 ### no self-spawning

### fisheries selectivity
sel. <- FLife::dnormal(age = ages_full + fish/nseasons, params = params)
#if (dims(sel.)["min"] == 0) sel.[1,,, 1] <- 0 ### recruits not fished
units(sel.) <- "f"

### natural mortality
### define Gislason natural mortality function
### FLife::gislason is wrong
gisM <- function(length, params) {
  exp(params["m1"] %+% (params["m2"] %*% log(length)) %+% 
        (params["m3"] %*% log(params["linf"])) %+% log(params["k"]))
}
### /seasons because M is annual value
m. <- gisM(length = lengths_mid, params = params)/range$seasons

### create FLBRP template
brp <- FLBRP(stock.wt        = weights_stk,
             landings.wt     = weights_catch,
             discards.wt     = weights_catch,
             bycatch.wt      = weights_catch,
             mat             = mat.,
             m               = m.,
             landings.sel    = sel.,
             discards.sel    = sel. %=% 0,
             bycatch.harvest = sel. %=% 0,
             harvest.spwn    = harvest.spwn,
             m.spwn          = m.spwn,
             availability    = weights_stk %=% 1,
             range           = unlist(range))

### set up unfished stock
# debugonce(N3)
### all metrics based on SSB=1000 when F=0
SSB0 <- 1000
Ri <- 1e+06
Ni <- N(object = flq, m = m., rec = Ri, eq_years = 100)
### SSB in 1st season (spawning time)
SSBi <- sum((Ni * mat. * weights_stk)[,,, 1]) 
spr0 <- SSBi/Ri
R0 <- SSB0 * 1/spr0

### unfished stock
N0 <- N(object = flq, m = m., rec = R0, eq_years = 100)
sum((N0 * mat. * weights_stk)[,,, 1])

### define recruitment model
model(brp) <- bevholt()$model
a <- (4 * params["s"] * R0)/(5 * params["s"] - 1)
b <- (SSB0 * (1 - params["s"]))/(5 * params["s"] - 1)
sr_annual <- as(brp, "FLSR")
params(sr_annual) <- FLPar(NA, dimnames = list(params = c("a","b"),
                                        iter = 1))
params(sr_annual)[,1] <- c(a, b)

### prepare OM stock
stk_annual <- FLStock(harvest = harvest(brp))
catch(stk_annual)[] <- catch(brp)
catch.n(stk_annual)[] <- catch.n(brp)
catch.wt(stk_annual)[] <- catch.wt(brp)
discards(stk_annual)[] <- discards(brp)
discards.n(stk_annual)[] <- discards.n(brp)
discards.wt(stk_annual)[] <- discards.wt(brp)
landings(stk_annual)[] <- landings(brp)
landings.n(stk_annual)[] <- landings.n(brp)
landings.wt(stk_annual)[] <- landings.wt(brp)
stock(stk_annual)[] <- stock(brp)
stock.n(stk_annual)[] <- stock.n(brp)
stock.wt(stk_annual)[] <- stock.wt(brp)
m(stk_annual)[] <- m(brp)
mat(stk_annual)[] <- mat(brp)
harvest(stk_annual)[] <- harvest(brp)
harvest.spwn(stk_annual)[] <- harvest.spwn(brp)
m.spwn(stk_annual)[] <- m.spwn(brp)
range(stk_annual)[] <- c(0, 5, 5, 1, 101, 1, 3)
#stk <- as(brp, "FLStock")[, 1:101]
for (y in 1) stock.n(stk_annual)[, y] <- N0 # do not fill later years -> weird things happen
stock(stk_annual) <- computeStock(stk_annual)
landings.n(stk_annual) <- discards.n(stk_annual) <- catch.n(stk_annual) <- 0
landings.n(stk_annual)[, 2:101] <- 1e-6
catch.n(stk_annual)[, 2:101] <- 1e-6
catch(stk_annual) <- computeCatch(stk_annual)
landings(stk_annual) <- computeLandings(stk_annual)
discards(stk_annual) <- computeDiscards(stk_annual)
harvest(stk_annual)[, 1] <- 0
for (y in 2:101) harvest(stk_annual)[, y] <- sel.

### zero fishing
#debugonce(fwd, signature = c("FLStock", "missing", "fwdControl"))
ctrl <- fwdControl(data.frame(year = rep(2:101),
                              season = 1,
                              quant = "fbar",
                              value = 0.0))
stk_annual_fwd0 <- fwd(stk_annual, control = ctrl, sr = sr_annual, effort_max = 1e+9)
### fish 0.4
ctrl <- fwdControl(data.frame(year = rep(2:101),
                              season = 1,
                              quant = "fbar",
                              value = 0.4/4))
stk_annual_fwd0.1 <- fwd(stk_annual, control = ctrl, sr = sr_annual)
### fish 1
ctrl <- fwdControl(data.frame(year = rep(2:101),
                              season = 1,
                              quant = "fbar",
                              value = 1))
stk_annual_fwd1 <- fwd(stk_annual, control = ctrl, sr = sr_annual)
### fish 4
ctrl <- fwdControl(data.frame(year = rep(2:101),
                              season = 1,
                              quant = "fbar",
                              value = 4))
stk_annual_fwd4 <- fwd(stk_annual, control = ctrl, sr = sr_annual, effort_max = 1e+9)
### fish 10
ctrl <- fwdControl(data.frame(year = rep(2:101),
                              season = 1,
                              quant = "fbar",
                              value = 10))
stk_annual_fwd_ <- fwd(stk_annual, control = ctrl, sr = sr_annual, effort_max = 1e+9)

plot(window(FLStocks("F=0" = stk_annual_fwd0, "F=0.4" = stk_annual_fwd0.1, 
                     "F=1" = stk_annual_fwd1,
                     "F=4" = stk_annual_fwd4,
                     "F=10" = stk_annual_fwd_), end = 10)) +
  theme_bw()


### find MSY - annual model
trace_env <- new.env()
F_annual(om_stk = stk_annual, om_sr = sr_annual, trace_env = trace_env, target = 0)
F_annual <- function(om_stk, om_sr, target, trace_env, seasons = 1) {
  res_trace_i <- mget("res_trace", envir = trace_env, 
                      ifnotfound = FALSE)$res_trace
  if (isFALSE(res_trace_i)) res_trace_i <- list()
  if (isTRUE(target %in% sapply(res_trace_i, function(x) x$target))) {
    res_list <- res_trace_i[[which(target == sapply(res_trace_i, 
                                                    function(x) x$target))[[1]]]]
    run <- FALSE
  } else {
    run <- TRUE
  }
  if (isTRUE(run)) {
    ctrl <- fwdControl(data.frame(year = rep(2:101, each = length(seasons)),
                                  season = seasons,
                                  quant = "fbar",
                                  value = target/length(seasons)))
    res <- fwd(om_stk, control = ctrl, sr = om_sr, effort_max = 1e+9)
    stat_yrs <- 91:100
    res_list <- list(
      target = target,
      TSB = median(tsb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      SSB = median(ssb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      Catch = median(apply(catch(res)[, ac(stat_yrs)], 2, sum), na.rm = TRUE),
      Fbar = median(apply(fbar(res)[, ac(stat_yrs)], 2, sum), na.rm = TRUE),
      Rec = median(rec(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      objective = median(apply(catch(res)[, ac(stat_yrs)], 2, sum), 
                                 na.rm = TRUE))
  }
  res_trace_i <- append(res_trace_i, list(res_list))
  res_trace_i <- unique(res_trace_i)
  assign(value = res_trace_i, x = "res_trace", envir = trace_env)
  print(c(unlist(res_list)))
  return(res_list$objective)
}

vals <- c(seq(0, 10, 0.5))
runs_annual <- foreach(target = vals) %do% {
  F_annual(om_stk = stk_annual, om_sr = sr_annual, trace_env = trace_env,
           target = target)
}
res_trace <- mget("res_trace", envir = trace_env)
runs_annual <- as.data.frame(do.call(bind_rows, res_trace))
MSY_annual <- optimise(f = F_annual, seasons = 1,
                   om_stk = stk_annual, om_sr = sr_annual, trace_env = trace_env,
                   interval = c(0, 10),
                   lower = 0, upper = 10,
                   maximum = TRUE,
                   tol = 0.00001)
res_trace <- mget("res_trace", envir = trace_env)
runs_annual <- as.data.frame(do.call(bind_rows, res_trace))


### find MSY - seasonal model (run again)
stk_seasonal <- readRDS("input/san/stk_seasonal_deterministic.rds")
sr_seasonal <- readRDS("input/san/sr_seasonal_deterministic.rds")
trace_env_seasonal <- new.env()
MSY_seasonal <- optimise(f = F_annual, 
                         om_stk = stk_seasonal, om_sr = sr_seasonal, seasons = 1:4,
                         trace_env = trace_env_seasonal,
                         interval = c(0, 10),
                         lower = 0, upper = 10,
                         maximum = TRUE,
                         tol = 0.00001)
# F_annual(om_stk = stk_seasonal, om_sr = sr_seasonal, seasons = 1:4,
#          trace_env = trace_env_seasonal, target = 0)
vals_seasonal <- c(seq(0, 10, 0.5))
runs_seasonal <- foreach(target = vals_seasonal) %do% {
  F_annual(om_stk = stk_seasonal, om_sr = sr_seasonal, seasons = 1:4,
           trace_env = trace_env_seasonal,
           target = target)
}
res_trace_seasonal <- mget("res_trace", envir = trace_env_seasonal)
runs_seasonal <- as.data.frame(do.call(bind_rows, res_trace_seasonal))


### annual MSY projection
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 1),
                              season = 1,
                              quant = "fbar",
                              value = MSY_annual$maximum/1))
MSY_run_annual <- fwd(stk_annual, control = ctrl, sr = sr_annual, effort_max = 1e+9)
plot(MSY_run_annual)
### seasonal MSY projection
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = MSY_seasonal$maximum/4))
MSY_run_seasonal <- fwd(stk_seasonal, control = ctrl, sr = sr_seasonal, effort_max = 1e+9)
plot(MSY_run_seasonal)
### zero F projections
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 1),
                              season = 1,
                              quant = "fbar",
                              value = 0))
zero_run_annual <- fwd(stk_annual, control = ctrl, sr = sr_annual, effort_max = 1e+9)
plot(zero_run_annual)
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 0))
zero_run_seasonal <- fwd(stk_seasonal, control = ctrl, sr = sr_seasonal, effort_max = 1e+9)
plot(zero_run_seasonal)


plot(FLStocks(annual = MSY_run_annual, seasonal = MSY_run_seasonal))
plot(window(FLStocks(annual = MSY_run_annual, seasonal = MSY_run_seasonal),
            end = 10))


### combine annual and seasonal projections
qnts <- FLQuants(MSY_annual_rec = rec(MSY_run_annual),
                 MSY_seasonal_rec = rec(MSY_run_seasonal),
                 MSY_annual_ssb = ssb(MSY_run_annual),
                 MSY_seasonal_ssb = ssb(MSY_run_seasonal),
                 MSY_annual_tsb = tsb(MSY_run_annual),
                 MSY_seasonal_tsb = tsb(MSY_run_seasonal),
                 MSY_annual_catch = catch(MSY_run_annual),
                 MSY_seasonal_catch = catch(MSY_run_seasonal),
                 MSY_annual_fbar = fbar(MSY_run_annual),
                 MSY_seasonal_fbar = fbar(MSY_run_seasonal),
                 zero_annual_rec = rec(zero_run_annual),
                 zero_seasonal_rec = rec(zero_run_seasonal),
                 zero_annual_ssb = ssb(zero_run_annual),
                 zero_seasonal_ssb = ssb(zero_run_seasonal),
                 zero_annual_tsb = tsb(zero_run_annual),
                 zero_seasonal_tsb = tsb(zero_run_seasonal),
                 zero_annual_catch = catch(zero_run_annual),
                 zero_seasonal_catch = catch(zero_run_seasonal),
                 zero_annual_fbar = fbar(zero_run_annual),
                 zero_seasonal_fbar = fbar(zero_run_seasonal)
                 )
qnts <- window(qnts, end = 12)
qnts_df <- as.data.frame(qnts)
qnts_df <- qnts_df %>%
  mutate(season = as.numeric(as.character(season)),
         year_season = year + (season - 1)/4) %>%
  separate(qname, into = c("target", "model", "quant"), sep = "_") %>%
  mutate(target = factor(target, levels = c("zero", "MSY"),
                        labels = c("Zero fishing", "MSY")),
         model = factor(model, levels = c("annual", "seasonal"),
                        labels = c("Annual model", "Seasonal model")),
         quant = factor(quant, levels = c("rec", "ssb", "tsb", "catch", "fbar"),
                        labels = c("Recruits", "SSB [t]", "TSB [t]", 
                                   "Catch [t]", "F (ages 1-3)"))) %>%
  mutate(data = ifelse(quant == "Recruits" & season != "1", 0, data))

qnts_df %>%
  ggplot(aes(x = year_season - 1, y = data, 
             colour = model, linetype = target, shape = target)) +
  geom_line(size = 0.2) +
  geom_point(size = 0.5) +
  #facet_wrap(~ quant, scales = "free_y", strip.position = "left", ncol = 3) +
  facet_grid(quant ~ model, scales = "free_y", switch = "y") +
  coord_cartesian(ylim = c(0, NA), xlim = c(0, 6)) +
  scale_x_continuous("year", breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_linetype_manual(name = "", values = c("solid", "2121")) +
  #scale_colour_discrete("") +
  scale_colour_brewer(name = "", palette = "Set1") +
  scale_shape("") +
  theme_bw(base_size = 8)  +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside",
        axis.title.y = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(1, "lines"), 
        legend.spacing.y = unit(0, "lines")
        #legend.position = c()
        )
ggsave("output/plots/san_OM_MSY_proj_annual_seasonal.png",
       width = 16, height = 12, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_OM_MSY_proj_annual_seasonal.pdf",
       width = 16, height = 12, units = "cm")
