\name{bootstrapUtil.checkForNegativeDataIfCV}
\alias{bootstrapUtil.checkForNegativeDataIfCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapUtil.checkForNegativeDataIfCV(group1, group2, Func)
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
function (group1, group2, Func) 
{
    if (Func == "cv") {
        if (any(c(group1, group2) < 0)) {
            stop("CV IS NONSENSICAL FOR DATA WITH NEGATIVE VALUES\n ...USE MEAN OR MEDIAN")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
