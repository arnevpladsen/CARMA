\name{CARMA}
\alias{CARMA}
\title{
Calculate CARMA scores
}
\description{
This is the main function in the \code{CARMA} package. The main input to the function is a list of segmented integer copy number profiles. The \code{CARMA} function calculates six CARMA scores within each defined genomic region.
}
\usage{
CARMA(segment.list, region.type = c("arm", "chrom", "bin"), bin.width = 5 * 10^7, overlap = T, hg.version = c("hg38","hg19"), gender = NULL, exclude.x.chrom = T, normalization.set = NULL)
}
\arguments{
  \item{segment.list}{
a list containing data frames with segmented integer copy number data. See \code{head(test.data[[1]])} for the format of the copy number data. The data frames shold contain 5 columns: chromosome, start position, end position, copynumber in major allele and copy number in minor allele.
}
  \item{region.type}{
the type of region in which CARMA scores are to be calculated: \code{"arm"}, \code{"crom"} or \code{"bin"} corresponding to chromosome arm, whole chromosome or bin. Default set to \code{"arm"}.
}
  \item{bin.width}{
the width of the bin region. Only used if \code{region.type="bin"}.
}
  \item{overlap}{
if \code{overlap=TRUE}, segments within a chromosome arm will overlap to fit within the chromosome arm. If \code{overlap=FALSE}, each bin will widen out or be trimmed a little bit relative the defined bin width to fit the chromosome arm. Only used if \code{region.type="bin"}.
}
  \item{hg.version}{
the build version: \code{"hg38"} or \code{"hg19"}. Default set to \code{"hg38"}.
}
  \item{gender}{
a character vector of length samples giving gender. If \code{gender="male"}, the functions calculating AMP and DEL will correct for an expected lower value for the X chromosome in segments covering the X chromosome. This is done by dividing the ploidy by 2 when calculating AMP and DEL in these segments. If \code{gender=NULL} female gender is assumed.
}
  \item{exclude.x.chrom}{
if \code{exclude.x.chrom=TRUE}, the X chromosome is excluded from the analysis.
}
  \item{normalization.set}{
if \code{normalization.set=NULL} the dataset currently being analyzed is used for normalization. Optionally a different set of CARMA scores can be used as a normalization set. The format of the data would then be the direct output of the \code{CARMA} function.
}
}

\details{
For more details about the algorithm see linked paper under \code{source}. 
}
\value{
A list containing six objects:
  \item{raw_regional_CARMA_scores}{A list containing raw regional CARMA scores.}
  \item{raw_genome_wide_CARMA_scores}{A list containing raw aggregated genome wide CARMA scores.}
  \item{normalized_regional_CARMA_scores}{A list containing normalized regional CARMA scores.}
  \item{normalized_genome_wide_CARMA_scores}{A list containing normalized aggregated genome wide CARMA scores.}
  \item{start_stop_per_region}{A data frame containing start and end positions of CARMA regions.}
  \item{parameters}{A list containing the input parameters used when running the CARMA function.}

}
\references{
%\insertRef{carma:preprint}{carma}
}
\examples{
# Format of input copy number profiles
head(test.data[[1]])

# Calculate CARMA scores per chromosome arm
carma.res <- CARMA(test.data, region.type="arm", hg="hg19")

# Calculate CARMA scores per whole chromosome
carma.res <- CARMA(test.data, region.type="chrom", hg="hg19")

# Calculate CARMA scores within 3*10^7 base pairs long regions
carma.res <- CARMA(test.data, region.type="bin", bin.width=3*10^7, hg="hg19")
}

\source{
http://biorxiv.org/lookup/doi/10.1101/769356
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


