\name{bootstrap.multiConfidenceIntervals}
\alias{bootstrap.multiConfidenceIntervals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.multiConfidenceIntervals(data, Func, groupNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{groupNames}{
%%     ~~Describe \code{groupNames} here~~
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
function (data, Func, groupNames = NULL) 
{
    if (is.null(groupNames)) {
        groupNames = colnames(data)
    }
    values = c()
    ints = matrix(nrow = ncol(data), ncol = 2, dimnames = list(groupNames, 
        c("lower", "upper")))
    for (gp in 1:ncol(data)) {
        cat("  ", groupNames[gp], "\n", sep = "")
        temp = bootstrap.confidenceInterval(data[, gp], plots = F, 
            Func = Func, verbose = F)
        values = c(values, temp$value)
        ints[gp, ] = temp$conf.int
    }
    return(list(ints = ints, values = values))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
