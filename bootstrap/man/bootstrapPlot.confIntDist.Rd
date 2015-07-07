\name{bootstrapPlot.confIntDist}
\alias{bootstrapPlot.confIntDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.confIntDist(pseudo, Func, trials, value, conf.int, plotParams, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pseudo}{
%%     ~~Describe \code{pseudo} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{conf.int}{
%%     ~~Describe \code{conf.int} here~~
}
  \item{plotParams}{
%%     ~~Describe \code{plotParams} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pseudo, Func, trials, value, conf.int, plotParams, 
    ...) 
{
    col = plotParams$col
    border = plotParams$border
    line.col = plotParams$line.col
    hist(pseudo, breaks = length(pseudo)/100, col = col, border = border, 
        main = paste("bootstrapped ", Func, "\n(n=", trials, 
            ")", sep = ""), xlab = "", ...)
    abline(v = value, col = line.col)
    abline(v = conf.int[1], col = line.col, lty = "dashed")
    abline(v = conf.int[2], col = line.col, lty = "dashed")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
