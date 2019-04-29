\name{poisson.blaker.limits}
\alias{poisson.blaker.limits}

\title{
  Blaker's Poisson confidence limits
}

\description{
  Fast and accurate calculation of Blaker's Poisson confidence limits.
}
\usage{
poisson.blaker.limits(x, level = 0.95, tol = 1e-10, ...)
}

\arguments{
  \item{x}{
   number of events.
}
  \item{level}{
   confidence level.
}
  \item{tol}{
   numerical tolerance.
}
  \item{...}{
   additional arguments to be passed to \code{poisson.blaker.lower.limit};
   in fact, just \code{maxiter} (see \code{\link{BlakerCI-internal}}). 
}

}
\details{
   Note that the Blaker's \code{(1 - alpha)} confidence interval 
   is the convex hull of the set \code{C} of those points 
   where the acceptability function (Blaker (2000)) exceeds 
   level \code{alpha}. 
   When \code{C} is not connected, the algorithm is, analogously to
   \code{binom.blaker.limits} (see its details),  
   immune from leaving out short intervals and making thus
   the confidence intervals over-liberal.
}

\value{
   Length 2 vector -- the lower and upper confidence limits.
}
\references{
%% ~put references to the literature/web site here ~
Blaker, H. (2000) Confidence curves and improved exact confidence 
   intervals for discrete distributions.
   Canadian Journal of Statistics 28: 783-798.\cr
   (Corrigenda: Canadian Journal of Statistics 29: 681.)

%Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
%   and Blaker's exact tests.
%   Biostatistics. 11:373-374.

%Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
%   Intervals for Discrete Data. 
%   R Journal 2(1): 53-58.

%Klaschka, J. (2010). BlakerCI: An algorithm and R package
%   for the Blaker's binomial confidence limits calculation.
%   Technical report No. 1099, Institute of Computer Science,
%   Academy of Sciences of the Czech Republic,
%   \url{http://hdl.handle.net/11104/0195722}.
%\url{http://www.library.sk/i2/content.csg.cls?ictx=cav&repo=crepo1&key=99673653505}

Lecoutre, B. & Poitevineau J. (2014). New results for computing
   Blaker's exact confidence interval limits for usual
   one-parameter discrete distributions.
   Communications in Statistics - Simulation and Computation, 
   \url{http://dx.doi.org/10.1080/03610918.2014.911900}.
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}
}

\note{
%%  ~~further notes~~
   Package \code{exactci} by M. P. Fay includes another algorithm
   that calculates Blaker's Poisson confidence limits
   (see user-level function \code{poisson.exact} and internal function 
   \code{exactpoissonCI}). 
   
   Lecoutre & Poitevineau (2014) designed another algorithm 
   for the calculation of the Blaker's confidence limits.
   It is closely analogous to that of \code{poisson.blaker.limits}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\tabular{ll}{
\code{exactci:poisson.exact} \tab
      One of the options yields Blaker's
      limits. %\cr
}
}

\examples{
poisson.blaker.limits(3) # [1] 0.8176914 8.5597971

}



%}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
