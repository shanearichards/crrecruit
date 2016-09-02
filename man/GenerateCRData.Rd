\name{GenerateCRData}
\alias{GenerateCRData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate capture-recapture data according to the recruitment model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a capture-recapture data set that is consistent with the recruitment model. The function returns a list that includes summarises of the data needed when fitting the model. This function is useful for testing the fitting approach.
}
\usage{
d <- GenerateCRData(T, E, N, t.bar, sigma, alpha.0, alpha.1 = 0, beta.0, beta.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
%%     ~~Describe \code{T} here~~
Vector of length \code{J} containing sampling days.
}
  \item{E}{
%%     ~~Describe \code{E} here~~
Vector of length \code{J} containing effort during each sampling day.
}
  \item{N}{
%%     ~~Describe \code{N} here~~
The total number of recruits that are availbale for capture during the survey.
}
  \item{t.bar}{
%%     ~~Describe \code{t.bar} here~~
The day that recruitment peaks.
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
Duration of the recruitment period: approximately 95\% of recruitment occurs over 4*sigma days.
}
  \item{alpha.0}{
%%     ~~Describe \code{alpha.0} here~~
Baseline catch rate.
}
  \item{alpha.1}{
%%     ~~Describe \code{alpha.1} here~~
Age-dependent change in catch rate.
}
  \item{beta.0}{
%%     ~~Describe \code{beta.0} here~~
Baseline mortality rate.
}
  \item{beta.1}{
%%     ~~Describe \code{beta.1} here~~
Age-dependent change in mortality rate.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Summary of survey (days of sampling and effort during each survey)
T <- c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28) # sampling days
E <- c(1, 1, 1,  1,  1,  1,  1,  1,  1,  1) # sampling hours per event

# True parameter values
N       <- 250  # total number of recruiting animals
t.bar   <- 3.0  # peak day of recruitment
sigma   <- 3.0  # spread in recruitment
alpha.0 <- 0.15 # baseline catchability
alpha.1 <- 0.0  # age-depependent catchability
beta.0  <- 0.03 # baseline mortality rate
beta.1  <- 0.10 # age-dependent mortality rate

# Generate some capture-recapture data
Simulated <- GenerateCRData(T, E, N, t.bar, sigma, alpha.0, alpha.1,
  beta.0, beta.1)

# Plot the true number of catchable animals for each survey event
plot(x = T, y = Simulated$N.catchable.obs,
  xlab = "Day", ylab = "Catchable animals")

# set up necessary global variables for llRecruit
y   <- Simulated$y
f   <- Simulated$f
l   <- Simulated$l
T.F <- Simulated$T.F
T.L <- Simulated$T.L

# Can the true parameters be recovered?
guess <- c(t.bar, sigma, alpha.0, alpha.1, beta.0, beta.1)
fit <- optim(par = guess, fn = llRecruit,
  lower = c(0, 1.0, 0.10, -0.1, 0.0025, 0.001),
  upper = c(6, 5.0, 0.35,  0.1, 0.1000, 0.190),
  method = "L-BFGS-B",
  control = list(fnscale = -1, trace = 1, maxit = 100, REPORT = 20))

fit$par # display maximum-likelihood parameter estimates

I <- dim(y)[1] # number of unique animals caught
I/(1-PrNotDetect(fit$par)) # estimate of population size (N)

# compare true and estimated values
plot(x = guess, y = fit$par, xlab = "True value", ylab = "Fitted value")
abline(a = 0, b = 1, col = "grey")

## The function is currently defined as
function (T, E, N, t.bar, sigma, alpha.0, alpha.1 = 0, beta.0,
    beta.1)
{
    last.sample.day <- max(T)
    J <- length(E)
    T.F <- min(T) - 14
    T.L <- max(T)
    Y.days <- last.sample.day - T.F + 1
    Y.t <- T.F:last.sample.day
    recruit.day <- seq(from = T.F, to = T.L, by = 1)
    u.True <- exp(-0.5 * ((recruit.day - t.bar)/sigma)^2)
    u.True <- u.True/sum(u.True)
    U <- cumsum(u.True)
    age <- 0:Y.days
    if (beta.1 != 0) {
        pr.alive <- exp(-(beta.0/beta.1) * (exp(beta.1 * age) -
            1))
    }
    else {
        pr.alive <- exp(-beta.0 * age)
    }
    pr.dead <- 1 - pr.alive
    r.1 <- runif(n = N, min = 0, max = 1)
    r.2 <- runif(n = N, min = 0, max = 1)
    day.recruited <- rep(0, N)
    day.died <- rep(0, N)
    longevity <- rep(0, N)
    for (i in 1:N) {
        day.recruited[i] <- sum(r.1[i] > U) + T.F
        longevity[i] <- sum(r.2[i] >= pr.dead) - 1
        day.died[i] <- day.recruited[i] + longevity[i]
    }
    Y <- matrix(data = 0, nrow = N, ncol = Y.days)
    for (j in 1:Y.days) {
        t <- T.F - 1 + j
        for (i in 1:N) {
            if ((t >= day.recruited[i]) & (t <= day.died[i])) {
                Y[i, j] <- 1
            }
        }
    }
    N.catchable <- rep(0, Y.days)
    for (j in 1:Y.days) {
        N.catchable[j] <- sum(Y[, j])
    }
    N.catchable.obs <- rep(0, Y.days)
    for (j in 1:Y.days) {
        if (Y.t[j] \%in\% T) {
            N.catchable.obs[j] <- sum(Y[, j])
        }
        else {
            N.catchable.obs[j] <- NA
        }
    }
    N.catchable.obs <- N.catchable.obs[!is.na(N.catchable.obs)]
    Y.obs <- NULL
    age.obs <- NULL
    for (j in 1:Y.days) {
        if (Y.t[j] \%in\% T) {
            Y.obs <- cbind(Y.obs, Y[, j])
            age.obs <- cbind(age.obs, Y.t[j] - day.recruited)
        }
    }
    y <- Y.obs
    for (j in 1:J) {
        for (i in 1:N) {
            if (Y.obs[i, j] == 1) {
                pr.capture <- 1 - exp(-E[j] * alpha.0 * exp(alpha.1 *
                  age.obs[i, j]))
                if (runif(1) <= pr.capture) {
                  y[i, j] <- 1
                }
                else {
                  y[i, j] <- 0
                }
            }
        }
    }
    y.tmp <- y
    y <- NULL
    for (i in 1:N) {
        if (sum(y.tmp[i, ]) > 0) {
            y <- rbind(y, y.tmp[i, ])
        }
    }
    I <- dim(y)[1]
    f <- rep(0, I)
    l <- rep(0, I)
    for (i in 1:I) {
        found.yet <- FALSE
        for (j in 1:J) {
            if (y[i, j] == 1) {
                if (!found.yet) {
                  f[i] <- j
                  found.yet <- TRUE
                }
                l[i] <- j
            }
        }
    }
    return(list(y = y, I = I, J = J, f = f, l = l, T.F = T.F,
        T.L = T.L, u.True = u.True, pr.alive = pr.alive, N.catchable = N.catchable,
        N.catchable.obs = N.catchable.obs))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
