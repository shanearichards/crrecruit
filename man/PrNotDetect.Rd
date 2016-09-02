\name{PrNotDetect}
\alias{PrNotDetect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Probability animals will not be detected during the survey
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function that returns the probability a randomly chosen animal will not be detected during the survey, given the survey design. Model parameters are provided in a vector. Survey information must be stored in specific-named global variables (see Details).
}
\usage{
PrNotDetect(params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
%%     ~~Describe \code{params} here~~
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


Suppose the survey involves the capture of animals over \code{J} sampling events. The data must be stored and summarised by the following R-variables:
\describe{
  \item{E}{A double vector of size \code{J} quantifying the effort applied during each sampling event.}
  \item{T}{An integer vector of size \code{J} providing the day of each sampling event.}
  \item{T.F}{An integer indicating the first possible day of recruitment.}
  \item{T.L}{An integer (> \code{T.F}) indicating the last possible day of recruitment.}
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A single value.
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
guess <- c(3.0, 3.0, 0.2, 0.0, 0.03, 0.1)
PrNotDetect(guess)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
