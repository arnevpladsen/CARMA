\name{plotSample}
\alias{plotSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot sample specific CARMA scores
}
\description{
This function creates sample specific plots of normalized regional CARMA scores. The plots are saved as png files.}
\usage{
plotSample(carma, samples=c(1), plot.dir="~/Desktop/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma}{
the output file from the \code{CARMA} function containing regional CARMA scores.
}
  \item{samples}{
a vector giving the samples to be plotted. The vector can contain either numeric values corresponding to the order of samples included in the segment list used as input for the \code{CARMA} function, or give the actual sample IDs corresponding to the names given in the segment list.
}
  \item{plot.dir}{
the directory in which the plots are to be printed.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
png files containg plots of sample specific CARMA scores.
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
plotSample(carma, samples=c(1,2,3), plot.dir="~/Desktop/")
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
