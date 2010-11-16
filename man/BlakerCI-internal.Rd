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
(\code{\link{binom.blaker.lower.limit}}) 
and a single acceptance value -- unmodified or \sQuote{unimodalized}
(\code{\link{binom.blaker.acc.single.p}}).
}
\usage{
binom.blaker.lower.limit(x, n, level, tol = 1e-10)
binom.blaker.acc.single.p(x, n, p, type = "orig", acc.tol = 1e-10, output = "acc")
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
      point (single value between 0 and 1) where to calculate the acceptance.
}
  \item{type}{
      "orig", or "unimod" -- either unmodified, or unimodalized acceptance
      (see \code{\link{binom.blaker.acc}}).      
}
  \item{acc.tol}{
      numerical tolerance (for the acceptance values when \code{type = "unimod"}). 
}
  \item{output}{
      the acceptance value output (\code{output = "acc"}, the default)
      or, instead, an auxiliary integer-valued parameter \code{q1}, used
      for testing whether points belong to the same continuous segment 
      of the acceptance function
      (\code{output = "q1"}), or both (\code{output = "both"}).
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
For \code{\link{binom.blaker.lower.limit}}, a single number -- 
the lower confidence limit. \cr
For \code{\link{binom.blaker.acc.single.p}} -- depending on the
\code{output} parameter -- a single acceptance value, or a single
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