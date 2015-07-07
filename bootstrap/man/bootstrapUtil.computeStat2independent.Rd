\name{bootstrapUtil.computeStat2independent}
\alias{bootstrapUtil.computeStat2independent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.computeStat2independent(group1, group2, Func = c("mean.diff", "median.diff", "cv.diff", "t.test", "wilcox.test", "ks.test"), verbose, dataDescriptor, replace)
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
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
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
function (group1, group2, Func = c("mean.diff", "median.diff", 
    "cv.diff", "t.test", "wilcox.test", "ks.test"), verbose, 
    dataDescriptor, replace) 
{
    choiceVec = c("mean.diff", "median.diff", "cv.diff _DO_NOT_USE_", 
        "t.test", "wilcox.test", "ks.test")
    title = "Please select one of the following statistics to use:"
    err = "SOMETHING IS WRONG IN bootstrapUtil.computeStat2independent(), GET AUSTIN"
    if (length(Func) > 1) {
        checkArg = menu(choices = choiceVec, title = title)
        Func = choiceVec[checkArg]
        cat("Using ", Func, " statistic...\n", sep = "")
    }
    else if (length(Func) == 1) {
        checkArg = pmatch(Func, choiceVec)
        if (is.na(checkArg)) {
            checkArg = menu(choices = choiceVec, title = title)
        }
        Func = choiceVec[checkArg]
    }
    else {
        stop(err)
    }
    if (checkArg \%in\% 1:3) {
        Func = gsub(".diff", "", Func)
        bootstrapUtil.checkForNegativeDataIfCV(group1, group2, 
            Func)
        stat = eval(call(Func, group1)) - eval(call(Func, group2))
        if (verbose) {
            cat("...........................................\n")
            cat("Testing: ", dataDescriptor, "\n\n", sep = "")
            cat("Building null distribution...\n")
            cat(" Test statistic: difference of group ", Func, 
                "s = ", stat, "\n", sep = "")
            cat(" Resampling with replacement: ", replace, "\n\n", 
                sep = "")
        }
    }
    else if (checkArg \%in\% 4:6) {
        Func = gsub(".test", "", Func)
        if (checkArg == 4) {
            stat = bootstrapUtil.tStat(group1, group2)
        }
        else if (checkArg == 5) {
            stat = bootstrapUtil.wilcoxStat(group1, group2)
        }
        else if (checkArg == 6) {
            stat = bootstrapUtil.ksStat(group1, group2)
        }
        if (verbose) {
            cat("...........................................\n")
            cat("Testing: ", dataDescriptor, "\n\n", sep = "")
            cat("Building null distribution...\n")
            cat(" Test statistic: ", Func, " = ", stat, "\n", 
                sep = "")
            cat(" Resampling with replacement: ", replace, "\n\n", 
                sep = "")
        }
    }
    else {
        stop(err)
    }
    return(list(stat = stat, Func = Func))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
