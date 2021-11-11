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
