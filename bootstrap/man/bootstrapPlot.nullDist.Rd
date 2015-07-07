\name{bootstrapPlot.nullDist}
\alias{bootstrapPlot.nullDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.nullDist(stat, abs.diffs, reflect, statsNULL, trials, Func, replace, col, border, col.line, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{
%%     ~~Describe \code{stat} here~~
}
  \item{abs.diffs}{
%%     ~~Describe \code{abs.diffs} here~~
}
  \item{reflect}{
%%     ~~Describe \code{reflect} here~~
}
  \item{statsNULL}{
%%     ~~Describe \code{statsNULL} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{border}{
%%     ~~Describe \code{border} here~~
}
  \item{col.line}{
%%     ~~Describe \code{col.line} here~~
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
function (stat, abs.diffs, reflect, statsNULL, trials, Func, 
    replace, col, border, col.line, ...) 
{
    hist_null = hist(statsNULL, breaks = floor(trials/100), plot = F)
    xmax = max(max(hist_null$breaks), abs(min(hist_null$breaks)))
    if (all(statsNULL > 0)) {
        xlim = c(0, xmax)
    }
    else if (all(statsNULL < 0)) {
        xlim = c(-xmax, 0)
    }
    else {
        xlim = c(-xmax, xmax)
    }
    if (Func \%in\% c("mean", "median", "cv")) {
        xlab = paste("statistic (", Func, ".diff, replace = ", 
            substr(replace, 1, 1), ")", sep = "")
    }
    else {
        xlab = paste("statistic (", Func, ", replace = ", substr(replace, 
            1, 1), ")", sep = "")
    }
    plot(hist_null, xlim = xlim, xlab = xlab, main = paste("bootstrap null distribution\n(n=", 
        trials, ")", sep = ""), col = col, border = border, ...)
    if (abs.diffs) {
        abline(v = stat, col = col.line)
    }
    else {
        abline(v = stat, col = col.line)
        abline(v = reflect, col = col.line, lty = "dashed")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
