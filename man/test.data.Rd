\name{test.data}
\alias{test.data}
\docType{data}
\title{
Test data - copy number profiles
}
\description{
A test data file containing segmented copy number profiles from 100 breast cancer samples, for test running the CARMA algorithm. The test data are assembled from real tumor profiles, by randomly combining chromosome-specific copy number profiles from multiple samples to construct new whole genome copy number profiles.
}
\usage{test.data}
\format{
  List of 100 data frames, each containing the following 5 columns:
  \describe{
    \item{\code{chr}}{a charachter giving the chromosome}
    \item{\code{startpos}}{a numeric vector giving the hg19 start position of the segment}
    \item{\code{endpos}}{a numeric vector giving the hg19 end position of the segment}
    \item{\code{nMajor}}{a numeric vector giving the major allele integer copy number state for the segment}
    \item{\code{nMinor}}{a numeric vector giving the minor allele integer copy number state for the segment}
  }
 
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
test.results = CARMA(test.data, region.type="arm", hg.version="hg19")
carma.summary = summaryCARMA(test.results)
}
\keyword{datasets}
