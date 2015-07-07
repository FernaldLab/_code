\name{bootstrap.2independent}
\alias{bootstrap.2independent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.2independent(group1, group2, Func = "mean", trials = 10000, replace = T, groupNames = c("group1", "group2"), plots = T, plotNullDist = T, savePlot = F, dataDescriptor = NULL, col = "grey", border = "darkgrey", col.line = "red", jitter = 0.15, pch = 21, verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
}
  \item{groupNames}{
%%     ~~Describe \code{groupNames} here~~
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
}
  \item{plotNullDist}{
%%     ~~Describe \code{plotNullDist} here~~
}
  \item{savePlot}{
%%     ~~Describe \code{savePlot} here~~
}
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
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
  \item{jitter}{
%%     ~~Describe \code{jitter} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
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
function (group1, group2, Func = "mean", trials = 10000, replace = T, 
    groupNames = c("group1", "group2"), plots = T, plotNullDist = T, 
    savePlot = F, dataDescriptor = NULL, col = "grey", border = "darkgrey", 
    col.line = "red", jitter = 0.15, pch = 21, verbose = T, ...) 
{
    allParam = as.list(sys.frame(sys.nframe()))
    tmp = bootstrapUtil.checkNumericAndNAs(group1, group2, paired = F)
    group1 = tmp$group1
    group2 = tmp$group2
    boxNULL = tmp$boxNULL
    tmp = bootstrapUtil.computeStat2independent(group1 = group1, 
        group2 = group2, Func = Func, replace = replace, verbose = verbose, 
        dataDescriptor = dataDescriptor)
    stat = tmp$stat
    Func = tmp$Func
    statsNULL = c()
    for (trial in 1:trials) {
        if (verbose) {
            bootstrapUtil.printRunNumber(trial = trial, trials = trials)
        }
        pseudo_group1 = sample(boxNULL, length(group1), replace = replace)
        pseudo_group2 = sample(boxNULL, length(group2), replace = replace)
        pseudo_stat_null = bootstrapUtil.computeStat2independent(group1 = pseudo_group1, 
            group2 = pseudo_group2, Func = Func, verbose = F, 
            dataDescriptor = dataDescriptor, replace = replace)$stat
        statsNULL = c(statsNULL, pseudo_stat_null)
    }
    if (any(is.na(statsNULL))) {
        stop("NA IN NULL DISTRIBUTION, CHECK YOUR DATA")
    }
    ptemp = bootstrapUtil.compute2sidePval(statsNULL = statsNULL, 
        stat = stat, trials = trials)
    if (plots) {
        lims = bootstrapUtil.match2histAxisLims(group1, group2)
        if (plotNullDist) {
            if (savePlot) {
                tmp = paste(deparse(match.call()$group1), "-", 
                  deparse(match.call()$group2), "_", Func, sep = "")
                jpeg(file = paste(tmp, ".jpg", sep = ""), width = 8, 
                  height = 8, units = "in", quality = 100, type = "quartz", 
                  res = 150)
            }
            par(mfrow = c(2, 2))
            bootstrapPlot.groupHist(group = group1, groupName = groupNames[1], 
                Func = Func, dataDescriptor = dataDescriptor, 
                lims = lims, col = col, border = border, col.line = col.line)
            bootstrapPlot.groupHist(group = group2, groupName = groupNames[2], 
                Func = Func, dataDescriptor = dataDescriptor, 
                lims = lims, col = col, border = border, col.line = col.line)
            bootstrapPlot.nullDist(stat = stat, abs.diffs = F, 
                reflect = ptemp$reflect, statsNULL = statsNULL, 
                trials = trials, Func = Func, replace = replace, 
                col = col, border = border, col.line = col.line)
            outerArgs = allParam[names(allParam) \%in\% names(formals(bootstrapPlot.2groupBoxplot))]
            do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, 
                p = ptemp$p, paired = F, ...))
        }
        else {
            if (savePlot) {
                tmp = paste(deparse(match.call()$group1), "-", 
                  deparse(match.call()$group2), "_", Func, sep = "")
                jpeg(file = paste(tmp, ".jpg", sep = ""), width = 12, 
                  height = 4, units = "in", quality = 100, type = "quartz", 
                  res = 150)
            }
            par(mfrow = c(1, 3))
            bootstrapPlot.groupHist(group = group1, groupName = groupNames[1], 
                Func = Func, dataDescriptor = dataDescriptor, 
                lims = lims, col = col, border = border, col.line = col.line)
            bootstrapPlot.groupHist(group = group2, groupName = groupNames[2], 
                Func = Func, dataDescriptor = dataDescriptor, 
                lims = lims, col = col, border = border, col.line = col.line)
            outerArgs = allParam[names(allParam) \%in\% names(formals(bootstrapPlot.2groupBoxplot))]
            do.call(bootstrapPlot.2groupBoxplot, c(outerArgs, 
                p = ptemp$p, paired = F, ...))
        }
        if (savePlot) {
            dev.off()
        }
    }
    if (verbose) {
        cat("\np = ", ptemp$p, "\n", sep = "")
    }
    output = list(stat = stat, stat.reflect = ptemp$reflect, 
        p = ptemp$p, p.left = ptemp$p_left, p.right = ptemp$p_right, 
        null.dist = statsNULL, null.dist.mean = ptemp$midNULL, 
        data = setNames(list(group1, group2), groupNames), param = list(call.param = match.call(), 
            all.param = allParam))
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
