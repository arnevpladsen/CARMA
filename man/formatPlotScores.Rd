\name{formatPlotScores}
\alias{formatPlotScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format CARMA scores for use in plotting functions
}
\description{
This function will not be called directly by the user. The function is used by the \code{plotSample} and \code{plotDataset} functions to reformat CARMA scores to an appropriate format.
}
\usage{
formatPlotScores(scores, reg.ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{
a table with CARMA scores obtained from the output file produced by the \code{CARMA} algorithm, found in the \code{[["normalized_regional_CARMA_scores"]]} list.
}
  \item{reg.ref}{
a table giving chromosome as well as start and stop positions for each CARMA region.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of data frames, each containing start and stop positions for each CARMA score.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
