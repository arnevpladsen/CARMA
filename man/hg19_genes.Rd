\name{hg19_genes}
\alias{hg19_genes}
\docType{data}
\title{
hg19 gene reference positions
}
\description{
Gene annotation file with genomic start and end positions for all genes. The gene annotation data is exctracted from the "Comprehensive gene annotation" GTF file for hg19, which is based on hg38, release 32 (GRCh38.p13) and mapped to hg19, obtained from the GENCODE website.
}
\usage{hg19_genes}
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
https://www.gencodegenes.org/human/release_32lift37.html
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
head(hg19_genes)
}
\keyword{datasets}
