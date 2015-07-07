\name{bootstrapUtil.compute2sidePval}
\alias{bootstrapUtil.compute2sidePval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.compute2sidePval(statsNULL, stat, trials)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statsNULL}{
%%     ~~Describe \code{statsNULL} here~~
}
  \item{stat}{
%%     ~~Describe \code{stat} here~~
}
  \item{trials}{
%%     ~~Describe \code{trials} here~~
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
function (statsNULL, stat, trials) 
{
    midNULL = mean(statsNULL)
    reflect = midNULL - (stat - midNULL)
    if (stat < midNULL) {
        p_left = sum(statsNULL < stat)/trials
        p_right = sum(statsNULL > reflect)/trials
        p = p_left + p_right
    }
    else if (stat > midNULL) {
        p_right = sum(statsNULL > stat)/trials
        p_left = sum(statsNULL < reflect)/trials
        p = p_right + p_left
    }
    else {
        stop("Either statistic==0 or midNULL or something else is wrong...\n")
    }
    return(list(midNULL = midNULL, reflect = reflect, p_left = p_left, 
        p_right = p_right, p = p))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
