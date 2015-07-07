\name{bootstrapPlot.2groupBoxplot}
\alias{bootstrapPlot.2groupBoxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapPlot.2groupBoxplot(group1, group2, groupNames, paired, p, dataDescriptor = NULL, col, border, pch, jitter, Func, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
%%     ~~Describe \code{group1} here~~
}
  \item{group2}{
%%     ~~Describe \code{group2} here~~
}
  \item{groupNames}{
%%     ~~Describe \code{groupNames} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{dataDescriptor}{
%%     ~~Describe \code{dataDescriptor} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{border}{
%%     ~~Describe \code{border} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{jitter}{
%%     ~~Describe \code{jitter} here~~
}
  \item{Func}{
%%     ~~Describe \code{Func} here~~
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
function (group1, group2, groupNames, paired, p, dataDescriptor = NULL, 
    col, border, pch, jitter, Func, ...) 
{
    groupNames[1] = paste(groupNames[1], " (n=", (length(group1) - 
        sum(is.na(group1))), ")", sep = "")
    groupNames[2] = paste(groupNames[2], " (n=", (length(group2) - 
        sum(is.na(group2))), ")", sep = "")
    toPlot = c(group1, group2)
    grp = c(rep(groupNames[1], length(group1)), rep(groupNames[2], 
        length(group2)))
    boxLineType = Func
    if (boxLineType == "mean") {
        medlty = "blank"
    }
    else {
        medlty = "solid"
    }
    toPaste = bootstrapUtil.getPvalToPaste(p)
    data.lab = bootstrapUtil.checkDataDescriptorToGetLabel(dataDescriptor = dataDescriptor)
    boxplot(toPlot ~ grp, ylab = data.lab, medlty = medlty, main = toPaste, 
        ...)
    if (boxLineType == "mean") {
        segments(0.6, mean(group1, na.rm = T), 1.4, mean(group1, 
            na.rm = T), lwd = 3, col = "black")
        segments(1.6, mean(group2, na.rm = T), 2.4, mean(group2, 
            na.rm = T), lwd = 3)
    }
    stripchart(toPlot ~ grp, vertical = T, add = T, method = "jitter", 
        jitter = jitter, pch = pch, bg = col, ...)
    if (paired) {
        data = data.frame(group1, group2)
        for (row in 1:nrow(data)) {
            segments(1, data[row, 1], 2, data[row, 2], col = border)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
