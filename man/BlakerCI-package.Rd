\name{BlakerCI-package}
\alias{BlakerCI-package}
\alias{BlakerCI}
\docType{package}
\title{
Blaker's binomial confidence limits
}
\description{
Fast and accurate calculation of Blaker's binomial
confidence limits.
}
\details{
\tabular{ll}{
Package: \tab BlakerCI\cr
Type: \tab Package\cr
Version: \tab 1.0-1\cr
Date: \tab 2010-11-26\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}

Maintainer: Jan Klaschka \email{klaschka@cs.cas.cz}
}

\references{
Blaker, H. (2000) Confidence curves and improved exact confidence 
   intervals for discrete distributions.
   Canadian Journal of Statistics 28: 783-798. \cr
   (Corrigenda: Canadian Journal of Statistics 29: 681.)

Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
   and Blaker's exact tests.
   Biostatistics. 11:373-374.

Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
   Intervals for Discrete Data. 
   R Journal 2(1): 53-58.

Klaschka, J. (2010). On calculation of Blaker's binomial 
   confidence limits. COMPSTAT presentation (extended version).\cr
   \url{www.cs.cas.cz/~klaschka/c10/417_ext.pdf}
}

\keyword{package}
\keyword{htest}

%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}

\examples{
binom.blaker.limits(3,10) # [1] 0.08726443 0.61941066
}
