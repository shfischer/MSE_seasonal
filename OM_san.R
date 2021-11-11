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
# remotes::install_github("flr/FLasher", ref = "fix_seasons",
#                         INSTALL_opts = "--no-multiarch")
# remotes::install_github("shfischer/FLasher", ref = "fix_seasons",
#                         INSTALL_opts = "--no-multiarch")
# remotes::install_github("shfischer/FLasher", ref = "fix_seasons",
#                         INSTALL_opts = "--no-multiarch")
source("funs.R")


### ------------------------------------------------------------------------ ###
### load life-history parameters ####
### ------------------------------------------------------------------------ ###
stocks <- read.csv(file = "input/stocks.csv")
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
ggsave("output/plots/san_length.png",
       width = 8.5, height = 6, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_length.pdf",
       width = 8.5, height = 6, units = "cm")



### ------------------------------------------------------------------------ ###
### test 4 seasons ####
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


### test stock projection
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


# debugonce(fwd, signature = c("FLStock", "missing", "fwdControl"))


### zero fishing
ctrl <- fwdControl(data.frame(year = rep(2:101, each = 4),
                              season = 1:4,
                              quant = "fbar",
                              value = 0.0))
stk_fwd0 <- fwd(stk, control = ctrl, sr = sr)
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
ggsave("output/plots/san_check_seasonal_F.png",
       width = 17, height = 10, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_check_seasonal_F.pdf",
       width = 17, height = 10, units = "cm")
plot(window(FLStocks("F=0" = stk_fwd0, "F=0.4" = stk_fwd0.1, "F=1" = stk_fwd1,
                     "F=4" = stk_fwd5), start = 95, end = 101)) +
  theme_bw()

### find MSY ####
fishF <- function(stk, sr, fbar, return_all = FALSE) {
  #browser()
  ctrl_tmp <- fwdControl(data.frame(year = rep(2:101, each = 4),
                                    season = 1:4,
                                    quant = "fbar",
                                    value = fbar/dim(stk)[4]))
  stk_fwd_tmp <- fwd(stk, control = ctrl_tmp, sr = sr, effort_max = 1e+9)
  # plot(stk_fwd_tmp)
  if (isTRUE(return_all)) {
    return(list(TSB = median(tsb(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE),
                SSB = median(ssb(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE),
                Catch = median(apply(catch(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                               na.rm = TRUE),
                Fbar = median(apply(fbar(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                              na.rm = TRUE),
                Rec = median(rec(stk_fwd_tmp)[, ac(92:101),, 1], na.rm = TRUE)))
  } else {
    return(median(apply(catch(stk_fwd_tmp)[, ac(92:101)], 2, sum), 
                  na.rm = TRUE))
  }
}
fishF(stk = stk, sr = sr, fbar = 0, return_all = TRUE)
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

res <- optimise(f = fishF, stk = stk, sr = sr, return_all = FALSE,
                interval = c(0, 5),
                lower = 0, upper = 4,
                maximum = TRUE,
                tol = 0.00001)
MSY_vals <- fishF(stk = stk, sr = sr, fbar = res$maximum, return_all = TRUE)
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
ggsave("output/plots/san_MSY.png",
       width = 17, height = 6, units = "cm", dpi = 300, type = "cairo")
ggsave("output/plots/san_MSY.pdf",
       width = 17, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### expand stock ####
### ------------------------------------------------------------------------ ###
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
saveRDS(sr_fwd, file = "input/san/sr.rds")

### ------------------------------------------------------------------------ ###
### create fishing histories ####
### ------------------------------------------------------------------------ ###

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
                  deviances = residuals(sr_fwd),
                  effort_max = 1e+9)
plot(window(stk_fwd_ow, end = 100))
plot(window(stk_fwd_ow, start = 90, end = 100))
saveRDS(stk_fwd_ow, file = "input/san/stk_ow.rds")

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
saveRDS(stk_fwd_rnd, file = "input/san/stk_rnd.rds")






