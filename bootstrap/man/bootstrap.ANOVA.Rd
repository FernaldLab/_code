\name{bootstrap.ANOVA}
\alias{bootstrap.ANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.ANOVA(data = "matrix or dataframe with groups in columns", Func = "mean", absDiffs = T, trials = 10000, replace = F, plots = T, dataDescriptor = NULL, groupNames = NULL, pch = 21, col = "grey", border = "darkgrey", line.col = "red", verbose = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{absDiffs}{
%%     ~~Describe \code{absDiffs} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
}
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{groupNames}{
%%     ~~Describe \code{groupNames} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{border}{
%%     ~~Describe \code{border} here~~
}
  \item{line.col}{
%%     ~~Describe \code{line.col} here~~
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
function (data = "matrix or dataframe with groups in columns", 
    Func = "mean", absDiffs = T, trials = 10000, replace = F, 
    plots = T, dataDescriptor = NULL, groupNames = NULL, pch = 21, 
    col = "grey", border = "darkgrey", line.col = "red", verbose = T, 
    ...) 
{
    calcF = bootstrapUtil.calculateF(data = data, Func = Func, 
        absDiffs = absDiffs)
    Fstat = calcF$FstatInfo$Fstat
    groupLengths = calcF$groupInfo$groupLengths
    nGroups = calcF$groupInfo$nGroups
    boxNULL = unlist(calcF$groupInfo$data)
    boxNULL = boxNULL[!(is.na(boxNULL))]
    statsNULL = c()
    statNAs = c()
    if (verbose) {
        cat("...........................................\n")
        cat("Testing: ", dataDescriptor, "\n\n", sep = "")
        cat("Building null distribution...\n")
        cat(" \"F\" statistic: based on group ", Func, "s = ", 
            Fstat, "\n", sep = "")
        cat(" Resampling with replacement: ", replace, "\n\n", 
            sep = "")
    }
    for (trial in 1:trials) {
        if (verbose) {
            bootstrapUtil.printRunNumber(trials = trials, trial = trial)
        }
        pseudoNULL = matrix(nrow = max(groupLengths), ncol = nGroups)
        for (gp in 1:ncol(pseudoNULL)) {
            if (groupLengths[gp] < nrow(pseudoNULL)) {
                numNA = nrow(pseudoNULL) - groupLengths[gp]
                pseudoNULL[1:numNA, gp] = NA
                pseudoNULL[(numNA + 1):nrow(pseudoNULL), gp] = sample(boxNULL, 
                  groupLengths[gp], replace = replace)
            }
            else {
                pseudoNULL[, gp] = sample(boxNULL, groupLengths[gp], 
                  replace = replace)
            }
        }
        pseudoF = bootstrapUtil.calculateF(data = pseudoNULL, 
            Func = Func, absDiffs = absDiffs)$FstatInfo$Fstat
        if (is.infinite(pseudoF) | is.nan(pseudoF) | is.na(pseudoF)) {
            statNAs = c(statNAs, trial)
            next
        }
        statsNULL = c(statsNULL, pseudoF)
    }
    if (length(statNAs) > 0) {
        warning("Check $statNAs for run #s where pseudoF was NA, consider setting absDiffs=T")
    }
    if (absDiffs) {
        p = sum(statsNULL > Fstat)/trials
        p_left = NULL
        p_right = NULL
        ptemp = list(p = p)
    }
    else {
        ptemp = bootstrapUtil.compute2sidePval(statsNULL = statsNULL, 
            stat = Fstat, trials = trials)
    }
    if (verbose) {
        cat("\nComputing 95\% confidence intervals...\n")
    }
    citmp = bootstrap.multiConfidenceIntervals(data = data, Func = Func, 
        groupNames = groupNames)
    ints = citmp$ints
    values = citmp$values
    if (plots) {
        par(mfrow = c(1, 2))
        hist(statsNULL, xlim = c(min(statsNULL), max(statsNULL, 
            Fstat)), xlab = "", main = paste("bootstrap null distribution\n(n=", 
            trials, ")", sep = ""), col = col, border = border, 
            ...)
        if (absDiffs) {
            abline(v = Fstat, col = line.col)
        }
        else {
            abline(v = Fstat, col = line.col)
            abline(v = ptemp$reflect, col = line.col, lty = "dashed")
        }
        if (is.null(groupNames)) {
            groupNames = colnames(data)
        }
        else if (!is.null(groupNames) & length(groupNames) == 
            ncol(data)) {
            groupNames = groupNames
        }
        boxplot(data, names = groupNames, main = paste("p = ", 
            ptemp$p, sep = ""), border = "lightgrey", medlty = "blank", 
            boxwex = 0.5, ...)
        stripchart(as.data.frame(data), vertical = T, add = T, 
            pch = pch, bg = col, ...)
        coMat = matrix(nrow = ncol(data), ncol = 2)
        coMat[1, ] = c(0.75, 1.25)
        for (gp in 2:nrow(coMat)) {
            coMat[gp, ] = coMat[gp - 1, ] + 1
        }
        for (gp in 1:length(values)) {
            segments(coMat[gp, 1], values[gp], coMat[gp, 2], 
                values[gp], lwd = 1, col = line.col)
            segments(coMat[gp, 1], ints[gp, 1], coMat[gp, 2], 
                ints[gp, 1], lwd = 1, lty = "dashed", col = line.col)
            segments(coMat[gp, 1], ints[gp, 2], coMat[gp, 2], 
                ints[gp, 2], lwd = 1, lty = "dashed", col = line.col)
        }
    }
    output = list(stat = Fstat, data = data, null.dist = statsNULL, 
        p = ptemp$p, p.left = ptemp$p_left, p.right = ptemp$p_right, 
        conf.ints = as.data.frame(ints), calcF = calcF, parameters = match.call(), 
        statNAs = statNAs)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
