\name{binRegions}
\alias{binRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Define CARMA bins
}
\description{
This function will not be called directly by the user. The function defines start and end positions of regional bins, in which CARMA scores later will be calculated. It is used by the \code{binSegments} function.
}
\usage{
binRegions(hg, bin.width, overlap)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hg}{
  	a data frame containing start and end positions of chromosome arms: \code{hg19} or \code{hg38}.
  	}
  \item{bin.width}{
  	the width of the bin region.
  	}
  \item{overlap}{if \code{overlap=TRUE}, segments within a chromosome arm will overlap to fit within the chromosome arm. If \code{overlap=FALSE}, each bin will widen out or be trimmed a little bit relative the defined bin width to fit the chromosome arm.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing start and end positions of bin regions.
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
binRegions(hg19, bin.width=3*10^7, overlap=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
