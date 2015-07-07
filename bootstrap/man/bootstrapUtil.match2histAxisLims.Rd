\name{bootstrapUtil.match2histAxisLims}
\alias{bootstrapUtil.match2histAxisLims}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.match2histAxisLims(group1, group2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
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
function (group1, group2) 
{
    hist_group1 = hist(group1, breaks = length(group1), plot = F)
    hist_group2 = hist(group2, breaks = length(group2), plot = F)
    ymax = max(c(hist_group1$counts, hist_group2$counts))
    xmin = min(c(hist_group1$breaks, hist_group2$breaks))
    xmax = max(c(hist_group1$breaks, hist_group2$breaks))
    return(list(ymax = ymax, xmin = xmin, xmax = xmax))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
