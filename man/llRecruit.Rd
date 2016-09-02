\name{llRecruit}
\alias{llRecruit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns the log-likelihood of the recruitment model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function that returns the log-likelihood of the model of recruitment, survival and catchability, given capture-recapture data. Model parameters are provided in the function's single argument, which is a vector. Capture-recapture data must be stored in specific-named global variables (see Details).
}
\usage{
llRecruit(params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{x} here~~
    A vector containing six parameter values (see Details).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The input \code{params} contains the following values:
\enumerate{
  \item [tbar]   The day when recruitment peaks.
  \item [sigma]  The spread in recruitment (days). Approximately 95\% of recruitments occur during days [tbar - 2*sigma, tbar + 2*sigma].
  \item [alpha0] Baseline catch rate (per unit effort).
  \item [alpha1] Age-dependence in catch rate.
  \item [beta0]  Base-line mortality rate (per day).
  \item [beta1]  Age-dependence in mortality rate.
}

Recruitment is described by a Guassian peak centred on day \code{tbar} with spread characterised by \code{sigma}. Per-capita catch rate is given by: \code{alpha0*exp(alpha1*age)}, where \code{age} is the days since recruitment. Per-capita mortality rate is given by: \code{beta0*exp(beta1*age)}.


Suppose the survey involves the capture of \code{I} animals over \code{J} sampling events. The data must be stored and summarised by the following R-variables:
\describe{
  \item{y}{An integer matrix containing the capture data (\code{I} rows and \code{J} columns). Rows indicate the animal and columns indicate the sampling event.}
  \item{f}{An integer vector of size \code{I} identifying the sampling event when the animal was first caught. Values lie in the range [1,\code{J}].}
  \item{l}{An integer vector of size \code{I} identifying the sampling event when the animal was last caught. Values lie in the range [1,\code{J}]. All element satisfy \code{l} >= \code{f}.}
  \item{E}{A double vector of size \code{J} quantifying the effort applied during each sampling event.}
  \item{T}{An integer vector of size \code{J} providing the day of each sampling event.}
  \item{T.F}{An integer indicating the first possible day of recruitment.}
  \item{T.L}{An integer (> \code{T.F}) indicating the last possible day of recruitment.}
}
Note that the vectors \code{f} and \code{l} are determined by the matrix \code{y} but must be calculated by the user.

This function is usually used in conjunction with the \code{optim} function (see Examples).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A double.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Shane A. Richards (shane.richards@csiro.au)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
guess <- c(3.0, 3.0, 0.2, 0.0, 0.03, 0.1) # dummy model parameters
ll <- llRecruit(guess) # calculate log-likelihood of the model

# find the maximum-likelihood parameter estimates
# box constraint method allows control over parameter estimates
#   age-dependent capture is turned off by restricting alpha1 (parameter 4)
# tracer flags are turned on to check for convergence
fit <- optim(par = guess, fn = llRecruit,
  lower = c(0, 1.0, 0.10, -0.001, 0.0025, 0.001),
  upper = c(6, 5.0, 0.35,  0.001, 0.1000, 0.190),
  method = "L-BFGS-B",
  control = list(fnscale = -1, trace = 1, maxit = 100, REPORT = 20))

fit$par # display maximum-likelihood parameter estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ecological statistics }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
