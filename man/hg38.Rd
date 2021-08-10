\name{hg38}
\alias{hg38}
\docType{data}
\title{
hg 38 chromosome arm start and end positions
}
\description{
A table with hg 38(GCA_000001405.15) start and end positions per chromosome arm for chromosome 1-22 and X. Data are exctracted from a hg38 cytoband position file, obtained from the University of California Santa Cruz, Genomics Institute website.
}
\usage{hg38}
\format{
  A data frame with the following 4 columns:  \describe{
    \item{\code{chr}}{a factor giving the chromosome}
    \item{\code{arm}}{a factor giving the chromosom arm}
    \item{\code{startpos}}{a numeric vector giving the hg38 start position}
    \item{\code{endpos}}{a numeric vector giving the hg 38 end position}
  }

}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
head(hg38)
}
\keyword{datasets}
