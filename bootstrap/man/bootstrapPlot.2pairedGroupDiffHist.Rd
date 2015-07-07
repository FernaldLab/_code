\name{bootstrapPlot.2pairedGroupDiffHist}
\alias{bootstrapPlot.2pairedGroupDiffHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.2pairedGroupDiffHist(dataDescriptor = NULL, diffs, Func, conditionNames, col, border, col.line, returnDataLabel = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{diffs}{
%%     ~~Describe \code{diffs} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{conditionNames}{
%%     ~~Describe \code{conditionNames} here~~
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
  \item{returnDataLabel}{
%%     ~~Describe \code{returnDataLabel} here~~
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
function (dataDescriptor = NULL, diffs, Func, conditionNames, 
    col, border, col.line, returnDataLabel = T, ...) 
{
    if (!is.null(dataDescriptor) & is.character(dataDescriptor)) {
        data.lab = dataDescriptor
    }
    else {
        data.lab = ""
        warning("Please label your data using arg 'dataDescriptor'\n  It's better for everyone")
    }
    hist.diffs = hist(diffs, breaks = length(diffs), plot = F)
    xmax = max(abs(hist.diffs$breaks))
    plot(hist.diffs, xlim = c(-xmax, xmax), main = paste("Individual differences with ", 
        Func, "\n(n=", length(diffs), ")", sep = ""), xlab = paste(data.lab, 
        " (", conditionNames[1], " - ", conditionNames[2], ")", 
        sep = ""), col = col, border = border, ...)
    abline(v = mean(diffs), col = col.line)
    if (returnDataLabel) {
        return(data.lab)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
