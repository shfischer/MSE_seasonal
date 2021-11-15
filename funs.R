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
### tmp ####
### ------------------------------------------------------------------------ ###

as_ <- function(from) {
        
        # GET 'quant' and dims
        
        qua <- quant(from[[1]])
        qdnms <- dimnames(from[[1]])[qua]
        itsq <- lapply(from, function(x) prod(dim(x[1,])[c(1,6)]))
        its <- max(unlist(itsq))
        
        # CONVERT to same quant
        from <- lapply(from, function(x) {
          if(dim(x)[1] == 1)
            dimnames(x) <- qdnms
          return(x)
        })
        
        # CONVERT
        df <- do.call("rbind", c(lapply(from, as.data.frame),
                                 make.row.names = FALSE))[,c(qua, "year", "iter", "data", "season")]
        df$qname <- rep(names(from), times=unlist(lapply(from, length)))
        
        # DEBUG as.data.frame(FLQuants) should accept qnames being equal if dims differ   
        # df <- as.data.frame(from)[,c(qua, "year", "iter",
        #   "data", "qname", "season")]
        
        # RESHAPE if min/max in quant
        if(any(df[,qua] %in% c("min", "max", "value"))) {
          df[,qua][df[,qua] == "all"] <- "value"
          df <- reshape(df, idvar = c("year", "iter", "qname", "season"),
                        timevar = qua, direction = "wide")
          names(df) <- gsub("data.", "", names(df))
          # or RENAME data as value
        } else {
          df[, qua] <- NULL
          names(df) <- sub("data", "value", names(df))
        }
        
        # RENAME qname to quant
        names(df) <- sub("qname", "quant", names(df))
        
        # DROP season if not used
        if(identical(unique(df$season), "all"))
          df$season <- NULL
        
        # NO ITERS
        if(its == 1) {
          
          target <- cbind(df[,-2], fishery=as.numeric(NA), catch=as.numeric(NA),
                          biol=1)
          
          return(fwdControl(target))
          
          # ITERS
        } else {
          
          target <- cbind(df[df$iter == df$iter[1],][,c('year', 'season', 'quant')])
          
          # ARRAY iters [targets, 3, iters]    
          iters <- array(NA, dim=c(dim(target)[1], 3, its),
                         dimnames=list(seq(dim(target)[1]), c("min", "value", "max"), 
                                       iter=seq(its)))
          
          # RESHAPE to assign from df
          # iters <- aperm(iters, c(3,1,2))
          iters[, "value", ] <- df$value
          if("min" %in% colnames(df))
            iters[, "min", ] <- df$min
          if("max" %in% colnames(df))
            iters[, "max", ] <- df$max
          # iters <- aperm(iters, c(2,3,1))
          
          # ADD fishery, catch and biol indices
          target <- cbind(target, fishery=as.numeric(NA), catch=as.numeric(NA),
                          biol=1)
          
          return(fwdControl(target=target, iters=iters))
        }
      }


### ------------------------------------------------------------------------ ###
### MSE loop with harvest rate ####
### ------------------------------------------------------------------------ ###
mse_hr <- function(om_stk, om_sr, 
                   yrs = 101:125, seasons = 1:4, 
                   hr = 0.1, lag = 0, interval = 1,
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
             by = interval)
  
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
                      dimnames = list(year = unique(tab$year[seq(id, id + interval - 1)]),
                                      season = unique(tab$season[seq(id, id + interval - 1)]),
                                      iter = dimnames(idx)$iter))
    advice[] <- rep(c(idx * hr), each = interval)
    ### set target
    ctrl <- as(FLQuants(catch = advice), "fwdControl")
    ### project OM
    om_stk <- fwd(om_stk, control = ctrl, sr = om_sr,
                  deviances = residuals(om_sr), effort_max = effort_max, 
                  maxF = maxF)
    #catch(om_stk)[, ac(yr),, ac(season)]
    #fbar(om_stk)[, ac(yr),, ac(season)]
    if (isTRUE(verbose)) cat("\n")
  }
  return(om_stk)
}


