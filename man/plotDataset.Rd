\name{plotDataset}
\alias{plotDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot summary of CARMA scores across the entire dataset or within subgroups in the dataset
}
\description{
This function plots regional CARMA scores across the whole dataset or within subgroups in the dataset. Each regional CARMA score is as a default presented as a median score across the dataset or within a group and can optionally be changed to reflect a percentile or a mean.
}
\usage{
plotDataset(carma, group=NULL, group.order=NULL, plot.titles=NULL, plot.path, n.plot.cols=1, max.color.score=0.7, FUN=median, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma}{
the output file from the \code{CARMA} function containing regional CARMA scores.
}
  \item{group}{
an optional argument giving a vector of length corresponding to the number of samples in the dataset. The vector assigns a label to each sample in order to plot each group separately. If \code{group=NULL} one summary plot of the entire dataset will be made. 
}
  \item{group.order}{
an optional argument in the form of a vector giving the order of each group to be plotted(from left to right and top to bottom), using the labels found in the \code{group} argument. Only relevant when the \code{group} argument is used.
}
  \item{plot.titles}{
an optional argument in the form of a vector giving the plot titles to be printed inside the circos plots. If \code{plot.titles=NULL} no plot titles will be printed.
}
  \item{plot.path}{
a file path for the png file. The file path must include a ".png" postfix.
}
  \item{n.plot.cols}{
number of circos plots per row.
}
  \item{max.color.score}{
a numeric value between 0 and 1 representing the dynamic range for the heat map. By choosing a low number, the color intensity for low scores is scaled up.
}
  \item{FUN}{
the function to be used when calculating CARMA scores across the dataset or within subgroups. \code{median} is chosen as a default, but can optionally be replaced by \code{mean} etc.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A png file containg a plot reflecting CARMA scores within the dataset.
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
carma <- CARMA(test.data, region.type="arm", hg.version="hg19")

#	Make plot of entire dataset
plotDataset(carma, plot.path="~/Desktop/carma_scores_entire_dataset.png")

#	Make plot of subgroups
plotDataset(carma, group=rep(c("pos", "neg"), 50), group.order=c("neg", "pos"), plot.titles = c("Negative", "Positive"), plot.path="~/Desktop/carma_scores_groups.png", n.plot.cols=1)

#	Make plots side by side by increasing argument "n.plot.cols"
plotDataset(carma, group=rep(c("pos", "neg"), 50), group.order=c("neg", "pos"), plot.titles = c("Negative", "Positive"), plot.path="~/Desktop/carma_scores_groups.png", n.plot.cols=3)

#	Intensify colors by decreasing the argument "max.color.score"
plotDataset(carma, group=rep(c("pos", "neg"), 50), group.order=c("neg", "pos"), plot.titles = c("Negative", "Positive"), plot.path="~/Desktop/carma_scores_groups.png", n.plot.cols=3, max.color.score=0.4)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
