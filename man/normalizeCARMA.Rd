\name{normalizeCARMA}
\alias{normalizeCARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize CARMA scores
}
\description{
This function will not be called directly by the user. The function is used by the \code{CARMA} function to calculate normalized regional CARMA scores. The function returns normalized CARMA scores per region for a list of samples, as well as normalized genome wide CARMA scores for a list of samples.
}
\usage{
normalizeCARMA(carma.res.reg, carma.res.genome, indx, region.type = "arm", bin.regions = NULL, exclude.x.chrom = T, normalization.set = NULL, quant = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma.res.reg}{
a list containing raw regional CARMA scores
}
  \item{carma.res.genome}{
a list containing aggregated regional CARMA scores (genome wide CARMA scores)
}
  \item{indx}{
names of the CARMA scores
}
  \item{region.type}{
type of region, either \code{"arm"}, \code{"chrom"} or \code{"bin"}.
}
  \item{bin.regions}{
if \code{region.type="bin"}, a data frame giving the start and end positions for the bins is required. It is generated using the \code{binRegions} function.
}
  \item{exclude.x.chrom}{
if \code{exclude.x.chrom=T}, the X chromosome is excluded from the analysis.
}
  \item{normalization.set}{
if \code{normalization.set=NULL} the dataset currently being analyzed is used for normalization. Optionally a different set of CARMA scores can be used as a normalization set. The format of the data would then be the direct output of the \code{CARMA} function.
}
  \item{quant}{
the percentile calculated within each region-specific CARMA score in the normalization set used for normalization.
}
}
\details{
Normalization is performed by first defining a normalization set, which by default is set to the dataset currently being analyzed. Optionally an external dataset can be provided for normalization. The input format of the normalization set would then be the direct output from the \code{CARMA} function. Next the 99th percentile for the CARMA score within each region is calculated in the normalization set. Each raw regional CARMA scores in each sample is then divided by the corresponding 99th percentile score.
}
\value{
Two lists of lengths corresponding to the number of samples. The first list (\code{[["regional"]]}) containes normalized regional CARMA scores. The second list (\code{[["genome"]]}) contains normalized aggregated genome wide CARMA scores.}
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
