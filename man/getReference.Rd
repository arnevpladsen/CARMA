\name{getReference}
\alias{getReference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Save CARMA input parameters
}
\description{
This function will not be called directly by the user. The function is used by the \code{CARMA} function and returns a table with genomic positions for the genomic region in which CARMA scores are calculated.
}
\usage{
getReference(region.type, hg, bin.width, overlap, exclude.x.chrom)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{region.type}{
the type of region chosen in which CARMA scores are calculated: \code{"arm"}, \code{"crom"} or \code{"bin"}. 
}
  \item{hg}{
a data frame containing start and end positions of chromosome arms: \code{hg19} or \code{hg38}.
}
  \item{bin.width}{
the width of the bin region.
}
  \item{overlap}{
if \code{overlap=TRUE}, segments within a chromosome arm will overlap to fit within the chromosome arm. If \code{overlap=FALSE}, each bin will widen out or be trimmed a little bit relative the defined bin width to fit the chromosome arm.
}
  \item{exclude.x.chrom}{
if \code{exclude.x.chrom=TRUE}, the X chromosome is excluded in the output file.
}
}
\details{

}
\value{
A data frame defining the start and end position for the genomic regions in which CARMA scores are calculated.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
