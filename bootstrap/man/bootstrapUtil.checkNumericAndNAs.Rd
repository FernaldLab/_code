\name{bootstrapUtil.checkNumericAndNAs}
\alias{bootstrapUtil.checkNumericAndNAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.checkNumericAndNAs(group1, group2, paired = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
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
function (group1, group2, paired = F) 
{
    not_numeric = !is.numeric(c(group1, group2))
    if (not_numeric) {
        stop("DATA CONTAINS NON-NUMERIC VALUES... FIX IT!\n")
    }
    if (paired) {
        nCheck = length(group1) != length(group2)
        if (nCheck) {
            stop("GROUPS ARE DIFFERENT SIZES...\n   DATA SHOULD BE PAIRED!!!")
        }
        NAcheck1 = is.na(group1)
        NAcheck2 = is.na(group2)
        removeMe = NAcheck1 | NAcheck2
        if (sum(removeMe) > 0) {
            group1 = group1[!removeMe]
            group2 = group2[!removeMe]
            warning("NAs in one/both conditions, corresponding data points removed from BOTH")
        }
        output = list(condition1 = group1, condition2 = group2)
    }
    else {
        NA_check = sum(is.na(c(group1, group2))) > 0
        if (NA_check) {
            group1 = group1[!is.na(group1)]
            group2 = group2[!is.na(group2)]
            boxNULL = c(group1, group2)
            warning("NAs removed from one or both groups, check your data")
        }
        else {
            boxNULL = c(group1, group2)
        }
        output = list(group1 = group1, group2 = group2, boxNULL = boxNULL)
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
