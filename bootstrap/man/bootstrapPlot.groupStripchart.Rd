\name{bootstrapPlot.groupStripchart}
\alias{bootstrapPlot.groupStripchart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.groupStripchart(data, interval, value, conf.int, Func, plotParams, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
}
  \item{value}{
%%     ~~Describe \code{value} here~~
}
  \item{conf.int}{
%%     ~~Describe \code{conf.int} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
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
function (data, interval, value, conf.int, Func, plotParams, 
    ...) 
{
    line.col = plotParams$line.col
    pch = plotParams$pch
    col = plotParams$col
    stripchart(data, vertical = T, pch = pch, bg = col, frame.plot = F, 
        main = paste(interval * 100, "\% conf.int around ", round(value, 
            1), "\n[", round(conf.int[1], 1), ", ", round(conf.int[2], 
            1), "]", sep = ""), ...)
    segments(0.8, value, 1.2, value, col = line.col)
    text(1.26, value, Func, ...)
    segments(0.9, conf.int[1], 1.1, conf.int[1], col = line.col, 
        lty = "dashed")
    segments(0.9, conf.int[2], 1.1, conf.int[2], col = line.col, 
        lty = "dashed")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
