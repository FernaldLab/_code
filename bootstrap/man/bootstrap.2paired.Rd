\name{bootstrap.2paired}
\alias{bootstrap.2paired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.2paired(condition1, condition2, trials = 10000, Func = "mean", plots = T, savePlot = F, abs.diffs = F, col = "grey", border = "darkgrey", col.line = "red", dataDescriptor = NULL, conditionNames = c("condition1", "condition2"), pch = 21, verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{condition1}{
%%     ~~Describe \code{condition1} here~~
}
  \item{condition2}{
%%     ~~Describe \code{condition2} here~~
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
  \item{savePlot}{
%%     ~~Describe \code{savePlot} here~~
}
  \item{abs.diffs}{
%%     ~~Describe \code{abs.diffs} here~~
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
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{conditionNames}{
%%     ~~Describe \code{conditionNames} here~~
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
function (condition1, condition2, trials = 10000, Func = "mean", 
    plots = T, savePlot = F, abs.diffs = F, col = "grey", border = "darkgrey", 
    col.line = "red", dataDescriptor = NULL, conditionNames = c("condition1", 
        "condition2"), pch = 21, verbose = T, ...) 
{
    allParam = as.list(sys.frame(sys.nframe()))
    tmp = bootstrapUtil.checkNumericAndNAs(condition1, condition2, 
        paired = T)
    condition1 = tmp$condition1
    condition2 = tmp$condition2
    diffs = condition1 - condition2
    stat = eval(call(Func, diffs))
    if (abs.diffs) {
        stat = abs(stat)
    }
    data = data.frame(condition1, condition2, diffs)
    names(data) = c(conditionNames[1], conditionNames[2], "diffs")
    if (verbose) {
        cat("...........................................\n")
        cat("Testing: ", dataDescriptor, "\n\n", sep = "")
        cat("Building null distribution...\n")
        cat(" Test statistic: ", Func, " of individual differences = ", 
            stat, "\n\n", sep = "")
    }
    statsNULL = c()
    for (trial in 1:trials) {
        if (verbose) {
            bootstrapUtil.printRunNumber(trials = trials, trial = trial)
        }
        signs = sample(c(-1, 1), length(diffs), replace = T)
        pseudo_diffs = signs * diffs
        pseudo_stat_null = eval(call(Func, pseudo_diffs))
        if (abs.diffs) {
            pseudo_stat_null = abs(pseudo_stat_null)
        }
        statsNULL = c(statsNULL, pseudo_stat_null)
    }
    if (any(is.na(statsNULL))) {
        stop("NA IN NULL DISTRIBUTION, CHECK YOUR DATA")
    }
    if (abs.diffs) {
        p = sum(statsNULL > stat)/trials
        ptemp = list(p = p)
    }
    else {
        ptemp = bootstrapUtil.compute2sidePval(statsNULL = statsNULL, 
            stat = stat, trials = trials)
    }
    if (plots) {
        par(mfrow = c(1, 3))
        outerArgs = allParam[names(allParam) \%in\% names(formals(bootstrapPlot.2pairedGroupDiffHist))]
        data.lab = do.call(bootstrapPlot.2pairedGroupDiffHist, 
            c(outerArgs, diffs = list(diffs), returnDataLabel = T, 
                ...))
        if (abs.diffs) {
            xlabNULL = paste("statistic: abs(", Func, " of individual differences)", 
                sep = "")
        }
        else {
            xlabNULL = paste("statistic: ", Func, " of individual differences", 
                sep = "")
        }
        hist(statsNULL, main = paste("bootstrap null distribution\n(n=", 
            trials, ")", sep = ""), xlab = xlabNULL, col = col, 
            border = border, breaks = floor(trials/100), ...)
        if (abs.diffs) {
            abline(v = stat, col = col.line)
        }
        else {
            abline(v = stat, col = col.line)
            abline(v = ptemp$reflect, col = col.line, lty = "dashed")
        }
        toPlot = data[, 1:2]
        if (ptemp$p == 0) {
            toPaste = paste("p < 1e-5 (n=", length(diffs), ")", 
                sep = "")
        }
        else {
            toPaste = paste("p = ", signif(ptemp$p, 2), " (n=", 
                length(diffs), ")", sep = "")
        }
        boxplot(toPlot, ylab = data.lab, main = toPaste, frame.plot = F, 
            ...)
        stripchart(toPlot, vertical = T, add = T, pch = pch, 
            bg = col, cex = 1.5)
        for (row in 1:nrow(data)) {
            segments(1, data[row, 1], 2, data[row, 2], col = border)
        }
    }
    if (abs.diffs) {
        output = list(stat = stat, p = ptemp$p, null.dist = statsNULL, 
            data = data, call.param = match.call(), all.param = allParam)
    }
    else {
        output = list(stat = stat, stat.reflect = ptemp$reflect, 
            p = ptemp$p, p.left = ptemp$p_left, p.right = ptemp$p_right, 
            null.dist = statsNULL, null.dist.mean = ptemp$midNULL, 
            data = data, call.param = match.call(), all.param = allParam)
    }
    if (verbose) {
        cat("\np = ", ptemp$p, "\n", sep = "")
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
