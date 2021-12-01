### recruitment model
N <- function(object, m, rec, eq_years = 100) {
  #browser()
  object[1, 1,, 1] <- rec
  for (a in seq(dim(object)[1])) {
    for (s in seq(dim(object)[4])) {
      ### recruitment 
      if (isTRUE(a == 1 & s == 1)) {
        next() 
        ### follow age through year
      } else if (isTRUE(s > 1 & a < dim(object)[1])) {
        object[a,,, s] <- object[a,,, s - 1] * exp(-m[a,,, s - 1])
        ### first season of new age
      } else if (isTRUE(s == 1 & a < dim(object)[1])) {
        object[a,,, s] <- object[a - 1,,, dim(object)[4]] * 
          exp(-m[a - 1,,, dim(object)[4]])
      ### last age (plusgroup)
      } else if (isTRUE(a == dim(object)[1])) {
        ### first season
        ### assume 100-year equilibrium
        if (isTRUE(s == 1)) {
          object[a,,, s] <- sum(c(object[a - 1,,, dim(object)[4]]) * 
            c(exp(-m[a - 1,,, dim(object)[4]]))^c(seq(eq_years)))
        } else {
          object[a,,, s] <- object[a,,, s - 1] * exp(-m[a,,, s - 1])
        }
        ### other seasons
      }
    }
  }
  return(object)
}

### ------------------------------------------------------------------------ ###
### fhist trajectories ####
### ------------------------------------------------------------------------ ###

fhist_random <- function(n_iter, yrs_hist, min = 0, max = 1) {
  
  start <- rep(0, n_iter)
  middle <- runif(n = n_iter, min = min, max = max)
  end <- runif(n = n_iter, min = min, max = max)
  df <- t(sapply(seq(n_iter), 
                 function(x) {
                   c(approx(x = c(1, yrs_hist/2), 
                            y = c(start[x], middle[x]), 
                            n = yrs_hist/2)$y,
                     approx(x = c(yrs_hist/2, yrs_hist + 1), 
                            y = c(middle[x], end[x]), 
                            n = (yrs_hist/2) + 1)$y[-1])
                 }))
  res <- FLQuant(c(t(df)), dimnames = list(year = seq(yrs_hist),
                                    iter = seq(n_iter)))
  return(res)
}

fhist_one_way <- function(n_iter, yrs_hist = 100, 
                          yrs_const = yrs_hist * 0.75, 
                          yrs_increase = yrs_hist * 0.25,
                          f0 = 0.5, fmax = 0.8) {
  
  ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
  fs <- rep(f0, yrs_const - 1)
  rate <- exp((log(fmax) - log(f0)) / (yrs_increase))
  fs <- c(fs, rate ^ (seq(yrs_increase)) * f0)
  res <- FLQuant(NA, dimnames = list(year = seq(yrs_hist),
                                     iter = seq(n_iter)))
  res[, 1] <- 0
  res[, -1] <- rep(fs, n_iter)
  return(res)
}

### ------------------------------------------------------------------------ ###
### MSE loop with harvest rate ####
### ------------------------------------------------------------------------ ###
mse_loop <- function(MP = "hr",
                     om_stk, om_sr, 
                     yrs = 101:125, seasons = 1:4, 
                     target = 0,
                     lag = 0, catch_interval = 1,
                     verbose = TRUE,
                     effort_max = 1e+9, maxF = 5) {
  #browser()
  tab <- expand.grid(year = as.numeric(dimnames(om_stk)$year),
                     season = as.numeric(dimnames(om_stk)$season))
  tab <- tab[order(tab$year, tab$season), ]
  row.names(tab) <- NULL
  tab$id <- seq(nrow(tab))
  
  ids <- seq(from = tab$id[tab$year == min(yrs) & tab$season == min(seasons)],
             to = tab$id[tab$year == max(yrs) & tab$season == max(seasons)], 
             by = catch_interval)
  
  for (id in ids) {
    yr <- tab$year[id]
    season <- tab$season[id]
    if (isTRUE(verbose)) 
      cat(paste0("year = ", tab$year[id], " - season = ", tab$season[id]))
    #yr = 101
    #season = 1
    #yr = yr+1
    #season = season+1
    ### get biomass at beginning of time step - requires forecast
    ctrl_tmp <- fwdControl(target = data.frame(quant = "fbar", year = yr,
                                               season = season, value = 0))
    om_stk_tmp <- fwd(om_stk, control = ctrl_tmp, sr = om_sr,
                      deviances = residuals(om_sr), effort_max = effort_max)
    ### get total biomass
    idx <- tsb(om_stk_tmp)[, ac(tab$year[id - lag]),, 
                           ac(tab$season[id - lag])]
    ### catch target
    advice <- FLQuant(NA,
                      dimnames = list(year = unique(tab$year[seq(id, id + catch_interval - 1)]),
                                      season = unique(tab$season[seq(id, id + catch_interval - 1)]),
                                      iter = dimnames(idx)$iter))
    if (identical(MP, "hr")) {
      advice[] <- rep(c(idx * target), each = catch_interval)
    } else if (identical(MP, "escapement")) {
      ### make sure catch is not negative
      ### if escapement biomass is below target, advice is zero 
      advice[] <- rep(pmax(c(idx - target), 0), each = catch_interval)
      ### fish excess throughout management period
      advice <- advice/catch_interval
    } else {
      stop("unknown MP")
    }
    ### set target
    ctrl <- as(FLQuants(catch = advice), "fwdControl")
    ### project OM
    om_stk_tmp <- fwd(om_stk, control = ctrl, sr = om_sr,
                      deviances = residuals(om_sr), effort_max = effort_max, 
                      maxF = maxF)
    ### insert only values from target years/seasons
    ### there is a bug in FLasher that previous years/seasons 
    ### might be overwritten with odd values...
    om_stk[, ac(ctrl@target$year),, ac(ctrl@target$season)] <- 
      om_stk_tmp[, ac(ctrl@target$year),, ac(ctrl@target$season)]
    #catch(om_stk)[, ac(yr),, ac(season)]
    #fbar(om_stk)[, ac(yr),, ac(season)]
    if (isTRUE(verbose)) cat("\n")
  }
  return(om_stk)
}

### ------------------------------------------------------------------------ ###
### optimise hr or escapement parameters ####
### ------------------------------------------------------------------------ ###

optimise_MP <- function(MP = "hr", om_stk = stk, om_sr = sr, 
                        yrs = 101:200, seasons = 1:4, 
                        target, lag = 0, catch_interval = 1,
                        verbose = FALSE, stat_yrs = 191:200,
                        return_all = FALSE, objective_stat = mean,
                        trace = FALSE, trace_env) {
  if (isTRUE(trace)) {
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
  } else {
    run <- TRUE
  }
  if (isTRUE(run)) {
    res <- mse_loop(MP = MP, om_stk = om_stk, om_sr = om_sr, 
                    yrs = yrs, seasons = seasons, 
                    target = target,
                    lag = lag, catch_interval = catch_interval,
                    verbose = verbose)
    res_list <- list(
      target = target,
      TSB = median(tsb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      SSB = median(ssb(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      Catch = median(apply(catch(res)[, ac(stat_yrs)], 2, sum), na.rm = TRUE),
      Fbar = median(apply(fbar(res)[, ac(stat_yrs)], 2, sum), na.rm = TRUE),
      Rec = median(rec(res)[, ac(stat_yrs),, 1], na.rm = TRUE),
      objective = objective_stat(apply(catch(res)[, ac(stat_yrs)], 2, sum), 
                                 na.rm = TRUE))
  }
  if (isTRUE(trace)) {
    res_add <- res_list
    res_trace_i <- append(res_trace_i, list(res_add))
    res_trace_i <- unique(res_trace_i)
    assign(value = res_trace_i, x = "res_trace", envir = trace_env)
  }
  print(c(unlist(res_list)))
  
  if (isTRUE(return_all)) {
    return(res_list)
  } else {
    return(res_list$objective)
    ### use mean so that if stock collapses in last year, objective value is
    ### reduced
  }
}

