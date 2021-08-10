\name{hg38_genes}
\alias{hg38_genes}
\docType{data}
\title{
hg38 gene reference positions
}
\description{
Gene annotation file with genomic start and end positions for all genes. The gene annotation data is exctracted from the "Comprehensive gene annotation" GTF file for hg38, release 32 (GRCh38.p13), obtained from the GENCODE website.
}
\usage{hg38_genes}
\format{
  This data frame contains the following 5 columns:
  \describe{
    \item{\code{chr}}{a factor giving the chromosome number}
    \item{\code{startpos}}{a numeric vector giving the hg38 start position of the gene}
    \item{\code{endpos}}{a numeric vector giving the hg38 end position of the gene}
    \item{\code{gene_id}}{a factor giving the Ensembl gene annotation ID}
    \item{\code{gene_name}}{a factor giving the gene name}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
https://www.gencodegenes.org/human/
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
head(hg38_genes)
}
\keyword{datasets}
