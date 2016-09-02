# a function that simulates multiple mark-recapture data
GenerateCRData <- function(T, E, N, t.bar, sigma, alpha.0, alpha.1 = 0.0, beta.0, beta.1) {
  last.sample.day <- max(T) # as a minimum set to the last observation day
  J <- length(E)            # sampling events
  T.F <- min(T) - 14        # earliest day of recruitment
  T.L <- max(T)             # latest day of recruitment
  Y.days <- last.sample.day - T.F + 1
  Y.t <- T.F:last.sample.day

  # calculate recruitment probability
  recruit.day <- seq(from = T.F, to = T.L, by = 1)
  u.True <- exp(-0.5*((recruit.day - t.bar)/sigma)^2)
  u.True <- u.True / sum(u.True)
  U <- cumsum(u.True) # cumulative recruitment distribution

  # calculate the cumulative mortality curve
  age <- 0:Y.days
  if (beta.1 != 0.0) {
    pr.alive <- exp(-(beta.0/beta.1)*(exp(beta.1*age)-1.0))
  } else {
    pr.alive <- exp(-beta.0*age)
  }
  pr.dead <- 1.0 - pr.alive

  #  generate demographic data
  r.1 <- runif(n = N, min = 0.0, max = 1.0) # random num for each individual
  r.2 <- runif(n = N, min = 0.0, max = 1.0) # random num for each individual
  day.recruited <- rep(0, N)
  day.died    <- rep(0, N)
  longevity   <- rep(0, N)
  for (i in 1:N) {
    day.recruited[i] <- sum(r.1[i] > U) + T.F # recruitment day
    longevity[i] <- sum(r.2[i] >= pr.dead) - 1 # whole days alive
    day.died[i] <- day.recruited[i] + longevity[i] # death day
  }

  # generate time series of unknown population size
  Y <- matrix(data = 0, nrow = N, ncol = Y.days) # alive matrix
  for (j in 1:Y.days) { # matrix time index
    t <- T.F - 1 + j # true day
    for (i in 1:N) {
      if ((t >= day.recruited[i]) & (t <= day.died[i])) {
        Y[i,j] <- 1
      }
    }
  }

  # calculate the number of catchable animals each day
  N.catchable <- rep(0,Y.days)
  for (j in 1:Y.days) {
    N.catchable[j] <- sum(Y[ ,j])
  }

  # filter sampling days
  N.catchable.obs <- rep(0,Y.days)
  for (j in 1:Y.days) {
    if (Y.t[j] %in% T) { # is this a survey day?
      N.catchable.obs[j] <- sum(Y[ ,j]) # yes, keep the data
    } else {
      N.catchable.obs[j] <- NA # no, remove the data for this day
    }
  }
  N.catchable.obs <- N.catchable.obs[!is.na(N.catchable.obs)] # rm NAs

  # generate time series of unknown population size for sampling days only
  Y.obs <- NULL
  age.obs <- NULL
  for (j in 1:Y.days) { # matrix time index
    if (Y.t[j] %in% T) {
      Y.obs <- cbind(Y.obs, Y[ ,j])
      age.obs <- cbind(age.obs, Y.t[j] - day.recruited)
    }
  }

  y <- Y.obs # observations
  for (j in 1:J) {
    for (i in 1:N) {
      if (Y.obs[i,j] == 1) {
        pr.capture <- 1.0 - exp(-E[j]*alpha.0*exp(alpha.1*age.obs[i,j]))
        if (runif(1) <= pr.capture) {
          y[i,j] <- 1
        } else {
          y[i,j] <- 0
        }
      }
    }
  }

  # remove non-captured animals from observations
  y.tmp <- y
  y <- NULL
  for (i in 1:N) {
    if (sum(y.tmp[i,]) > 0) {
      y <- rbind(y, y.tmp[i, ])
    }
  }

  I <- dim(y)[1] # number of animals marked

  # generate first and last sampling event vectors needed for fitting
  f <- rep(0, I)
  l <- rep(0, I)
  for (i in 1:I) {
    found.yet <- FALSE
    for (j in 1:J) {
      if (y[i,j] == 1) {
        if (!found.yet) {
          f[i] <- j
          found.yet <- TRUE
        }
        l[i] <- j
      }
    }
  }

  return (list(y = y, I = I, J = J, f = f, l = l, T.F = T.F, T.L = T.L,
    u.True = u.True, pr.alive = pr.alive,
    N.catchable = N.catchable, N.catchable.obs = N.catchable.obs))
}
