\name{bootstrapPlot.groupHist}
\alias{bootstrapPlot.groupHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.groupHist(group, groupName, Func, dataDescriptor = NULL, lims = NULL, col, border, col.line, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{groupName}{
%%     ~~Describe \code{groupName} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
}
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{lims}{
%%     ~~Describe \code{lims} here~~
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
function (group, groupName, Func, dataDescriptor = NULL, lims = NULL, 
    col, border, col.line, ...) 
{
    groupName = paste(groupName, " (n=", length(group), ")", 
        sep = "")
    hist_group = hist(group, breaks = length(group), plot = F)
    if (is.null(lims)) {
        ymax = max(hist_group$counts)
        xmin = min(hist_group$breaks)
        xmax = max(hist_group$breaks)
    }
    else if (!is.null(lims)) {
        if (is.list(lims) & is.numeric(unlist(lims))) {
            ymax = lims$ymax
            xmin = lims$xmin
            xmax = lims$xmax
        }
    }
    data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor)
    plot(hist_group, main = groupName, ylim = c(0, ymax), xlim = c(xmin, 
        xmax), col = col, border = border, xlab = data.lab, ...)
    if (Func \%in\% c("mean", "t")) {
        v = mean(group)
    }
    else if (Func \%in\% c("median", "cv", "wilcox", "ks")) {
        v = median(group)
    }
    abline(v = v, col = col.line)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
