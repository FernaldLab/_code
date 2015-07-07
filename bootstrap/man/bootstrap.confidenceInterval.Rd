\name{bootstrap.confidenceInterval}
\alias{bootstrap.confidenceInterval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.confidenceInterval(data, interval = 0.95, trials = 10000, Func = "mean", plots = T, showDist = T, savePlot = F, plotParams = list(col = "grey", border = "darkgrey", line.col = "red", pch = 21), verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
}
  \item{showDist}{
%%     ~~Describe \code{showDist} here~~
}
  \item{savePlot}{
%%     ~~Describe \code{savePlot} here~~
}
  \item{plotParams}{
%%     ~~Describe \code{plotParams} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (data, interval = 0.95, trials = 10000, Func = "mean", 
    plots = T, showDist = T, savePlot = F, plotParams = list(col = "grey", 
        border = "darkgrey", line.col = "red", pch = 21), verbose = T, 
    ...) 
{
    if (!is.numeric(data) | !is.vector(data)) {
        stop("CHECK THAT DATA IS A NUMERIC VECTOR")
    }
    if (any(is.na(data))) {
        data0 = data
        data = data[!(is.na(data))]
        warning(paste(sum(is.na(data0)), " NAs detected in data", 
            sep = ""))
    }
    value = eval(call(Func, data))
    if (verbose) {
        cat("....................................................................\n")
        cat("Computing ", interval * 100, "\% confidence interval around the ", 
            Func, "...\n", sep = "")
        cat(" Resampling...\n")
    }
    pseudo = c()
    for (trial in 1:trials) {
        if (verbose) {
            bootstrapUtil.printRunNumber(trials = trials, trial = trial)
        }
        pseudoData = sample(data, length(data), replace = T)
        pseudo = c(pseudo, eval(call(Func, pseudoData)))
    }
    half = (1 - interval)/2
    int = c(half, interval + half)
    conf.int = quantile(pseudo, int)
    if (plots) {
        if (showDist) {
            if (savePlot) {
                jpeg(file = paste(deparse(match.call()$data), 
                  ".jpg", sep = ""), width = 8, height = 4, units = "in", 
                  quality = 100, type = "quartz", res = 150)
            }
            par(mfrow = c(1, 2))
            bootstrapPlot.confIntDist(pseudo = pseudo, Func = Func, 
                trials = trials, value = value, conf.int = conf.int, 
                plotParams = plotParams, ...)
            bootstrapPlot.groupStripchart(data = data, interval = interval, 
                value = value, conf.int = conf.int, Func = Func, 
                plotParams = plotParams, ...)
        }
        else {
            if (savePlot) {
                jpeg(file = paste(deparse(match.call()$data), 
                  ".jpg", sep = ""), width = 4, height = 4, units = "in", 
                  quality = 100, type = "quartz", res = 150)
            }
            bootstrapPlot.groupStripchart(data = data, interval = interval, 
                value = value, conf.int = conf.int, Func = Func, 
                plotParams = plotParams, ...)
        }
        if (savePlot) {
            dev.off()
        }
    }
    output = list(value = value, conf.int = conf.int, dist = pseudo, 
        call.param = match.call())
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
