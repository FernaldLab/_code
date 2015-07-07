\name{bootstrapUtil.calculateF}
\alias{bootstrapUtil.calculateF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.calculateF(data, Func = "mean", absDiffs = T)
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
function (data, Func = "mean", absDiffs = T) 
{
    if (!is.data.frame(data) & !is.matrix(data)) {
        stop("DATA MUST BE IN MATRIX OR DATAFRAME")
    }
    if (mode(data) != "numeric") {
        stop("DATA MUST BE NUMERIC")
    }
    if (nrow(data) < ncol(data)) {
        warning("More groups than data points per group...\n   ARE YOU SURE?")
    }
    data = as.data.frame(data)
    nGroups = ncol(data)
    groupNums = apply(data, 2, get(Func), na.rm = T)
    groupNAs = apply(apply(data, 2, is.na), 2, sum)
    groupLengths = nrow(data) - groupNAs
    groupInfo = list(data = data, nGroups = nGroups, groupNums = groupNums, 
        groupNAs = groupNAs, groupLengths = groupLengths)
    grandTop = sum(groupLengths * groupNums)
    grandBot = sum(groupLengths)
    grandNum = grandTop/grandBot
    grandNumInfo = list(grandTop = grandTop, grandBot = grandBot, 
        grandNum = grandNum)
    diffs = grandNum - groupNums
    if (absDiffs) {
        diffs = abs(diffs)
    }
    numerator = sum(groupLengths * diffs)
    Fnum = list(diffs = diffs, numerator = numerator)
    toSum = list()
    sums = c()
    for (gp in 1:nGroups) {
        temp = data[, gp] - groupNums[gp]
        temp = temp[!is.na(temp)]
        if (absDiffs) {
            temp = abs(temp)
        }
        toSum[[gp]] = temp
        sums[gp] = sum(toSum[[gp]])
    }
    denominator = sum(sums)
    Fden = list(toSum = toSum, sums = sums, denominator = denominator)
    Fstat = numerator/denominator
    FstatInfo = list(Fnum = Fnum, Fden = Fden, Fstat = Fstat)
    out = list(groupInfo = groupInfo, grandNumInfo = grandNumInfo, 
        FstatInfo = FstatInfo, parameters = match.call())
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
