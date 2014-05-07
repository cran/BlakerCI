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
Version: \tab 1.0-4\cr
Date: \tab 2014-05-07\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}

\author{
Jan Klaschka \email{klaschka@cs.cas.cz}

Maintainer: Jan Klaschka \email{klaschka@cs.cas.cz}
}





%\references{
%Blaker, H. (2000) Confidence curves and improved exact confidence 
%   intervals for discrete distributions.
%   Canadian Journal of Statistics 28: 783-798. \cr
%   (Corrigenda: Canadian Journal of Statistics 29: 681.)
%
%Fay, M. P. (2010). Confidence intervals that Match Fisher's exact 
%   and Blaker's exact tests.
%   Biostatistics. 11:373-374.
%
%Fay, M.P. (2010a). Two-sided Exact Tests and Matching Confidence 
%   Intervals for Discrete Data. 
%   R Journal 2(1): 53-58.
%
%Klaschka, J. (2010). BlakerCI: An algorithm and R package
%   for the Blaker's binomial confidence limits calculation.
%   Technical report No. 1099, Institute of Computer Science,
%   Academy of Sciences of the Czech Republic.\cr
%\url{http://www.library.sk/i2/content.csg.cls?ictx=cav&repo=crepo1&key=99673653505}
%}


\section{Acknowledgement}{
This work was accomplished with institutional support RVO:67985807. 

Previous versions were supported by grant 205/09/1079 from 
The Czech Science Foundation.
%
%I am grateful to Petr Savicky for a valuable e-mail
%discussion of January 2010.
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
