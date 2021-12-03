### ------------------------------------------------------------------------ ###
### some seasonal projections ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
if (exists(x = "args_local")) args <- append(args, args_local)
print("arguments passed on to this script:")
print(args)
for (i in seq_along(args)) eval(parse(text = args[[i]]))

### ------------------------------------------------------------------------ ###
### set-up ####
### ------------------------------------------------------------------------ ###
### manual overwrite of library path for parallel workers...
.libPaths(c("renv/library/R-4.1/x86_64-conda-linux-gnu/", .libPaths()))
source("funs.R")
library(FLasher)
library(doParallel)
library(tidyverse)
library(cowplot)

# n_workers <- 10
cl <- makeCluster(n_workers)
registerDoParallel(cl)
. <- foreach(seq(length(cl))) %dopar% {
  .libPaths(c("renv/library/R-4.1/x86_64-conda-linux-gnu/", .libPaths()))
  source("funs.R")
  library(FLasher)
}

### ------------------------------------------------------------------------ ###
### HR with annual TAC - lags ####
### ------------------------------------------------------------------------ ###

if (identical(MP, "hr")) {

  stk_rnd_500 <- readRDS("input/san/stk_rnd_500.rds")
  sr_500 <- readRDS("input/san/sr_500.rds")
  hr_res_MSY_annual <- readRDS("input/san/hr_res_MSY_annual.rds")
  hr_target <- hr_res_MSY_annual$maximum ### annual MSY HR
  refpts <- readRDS("input/san/refpts.rds")
  
  iters_list <- split(seq(500), sort(seq(500) %% n_blocks))
  input_list <- lapply(iters_list, function(x) {
    list(stk = iter(stk_rnd_500, x), sr = iter(sr_500, x))
  })
  
  ### run for lag 0, 1, ... 8
  stk_lags <- foreach(lag_i = lag) %do% {
    print(paste0("lag=", lag_i))
    stk_junks <- foreach(input = input_list) %dopar% {
      # browser()
      res_i <- mse_loop(MP = "hr", force_seasonal = TRUE,
                        om_stk = input$stk, om_sr = input$sr, 
                        yrs = 101:200, seasons = 1:4, 
                        target = hr_target, lag = lag_i, catch_interval = 4,
                        verbose = FALSE)
      return(res_i)
    }
    stk_out <- stk_rnd_500
    for (i in seq_along(stk_junks)) 
      iter(stk_out, iters_list[[i]]) <- stk_junks[[i]]
    saveRDS(stk_out, paste0("output/san/seasonal/rnd_MSY_hr_annual_lag", lag_i,
                            ".rds"))
    print(paste0("finished lag=", lag_i))
  }
  
}

### ------------------------------------------------------------------------ ###
### Escapement strategy with annual TAC - lags ####
### ------------------------------------------------------------------------ ###

if (identical(MP, "hr")) {

  stk_rnd_500 <- readRDS("input/san/stk_rnd_500.rds")
  sr_500 <- readRDS("input/san/sr_500.rds")
  esc_res_MSY_annual <- readRDS("input/san/esc_res_annual_MSY.rds")
  esc_target <- esc_res_MSY_annual$maximum ### annual MSY HR
  refpts <- readRDS("input/san/refpts.rds")
  
  iters_list <- split(seq(500), sort(seq(500) %% n_blocks))
  input_list <- lapply(iters_list, function(x) {
    list(stk = iter(stk_rnd_500, x), sr = iter(sr_500, x))
  })
  
  ### run for lag 0, 1, ... 8
  stk_esc_lags <- foreach(lag_i = lag) %do% {
    print(paste0("starting lag=", lag_i))
    stk_junks <- foreach(input = input_list) %dopar% {
      # browser()
      res_i <- mse_loop(MP = "escapement", om_stk = input$stk, om_sr = input$sr, 
                        yrs = 101:200, seasons = 1:4, 
                        target = esc_target, lag = lag_i, catch_interval = 4,
                        verbose = FALSE)
      return(res_i)
    }
    stk_out <- stk_rnd_500
    for (i in seq_along(stk_junks)) 
      iter(stk_out, iters_list[[i]]) <- stk_junks[[i]]
    saveRDS(stk_out, paste0("output/san/seasonal/rnd_MSY_esc_annual_lag", lag_i,
                            ".rds"))
    print(paste0("finished lag=", lag_i))
  }
  
}



