\name{BlakerCI-internal}
\alias{BlakerCI-internal}
\alias{binom.blaker.lower.limit}
\alias{binom.blaker.acc.single.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions, not expected to be called by the user.
}
\description{
Calculation of the lower Blaker's confidence limit 
(\code{binom.blaker.lower.limit}) 
and a single acceptability value -- unmodified or \sQuote{unimodalized}
(\code{binom.blaker.acc.single.p}).
}
\usage{
binom.blaker.lower.limit(x, n, level, tol = 1e-10, maxiter=100)
binom.blaker.acc.single.p(x, n, p, type = "orig", acc.tol = 1e-10,
      output = "acc", maxiter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      number of successes.
}
  \item{n}{
      number of trials.
}
  \item{level}{
      confidence level.
}
  \item{tol}{
      numerical tolerance (for the confidence limit).
}
  \item{p}{
      point (single value between 0 and 1) where to calculate the acceptability.
}
  \item{type}{
      "orig", or "unimod" -- either unmodified, or unimodalized acceptability \cr
      (see \code{\link{binom.blaker.acc}}).      
}
  \item{acc.tol}{
      numerical tolerance (for the acceptability values when \code{type = "unimod"}). 
}
  \item{output}{
      the acceptability value output (\code{output = "acc"}, the default)
      or, instead, an auxiliary integer-valued parameter \code{q1}, used
      for testing whether points belong to the same continuous segment 
      of the acceptability function
      (\code{output = "q1"}), or both (\code{output = "both"}).
}
  \item{maxiter}{
      Maximum number of interval halving iterations during the search
      for a confidence limit (\code{binom.blaker.lower.limit}) or
      a discontinuity point of the acceptability function 
      (\code{binom.blaker.acc.single.p} with \code{type} \code{=} \code{"unimod"}). 
      When the required accuracy is not reached in \code{maxiter} steps -- 
      typically when too small \code{tol} or \code{acc.tol}
      exceeds capabilities of machine arithmetic -- last step's result
      is returned with warning.
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
For \code{binom.blaker.lower.limit}, a single number -- 
the lower confidence limit. \cr
For \code{binom.blaker.acc.single.p} -- depending on the
\code{output} parameter -- a single acceptability value, or a single
auxiliary integer, or both.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{binom.blaker.limits}}, \code{\link{binom.blaker.acc}}.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
