\name{genesCARMA}
\alias{genesCARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for mapping of CARMA scores to genes
}
\description{
This function assigns CARMA scores to genes.
}
\usage{
genesCARMA(carma, data.type = "normalized")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma}{
a list containing CARMA scores for each sample, obtained by running the \code{CARMA} function.
}
  \item{data.type}{
\code{"normalized"} or \code{"raw"}, corresponding to normalized or raw CARMA scores to be mapped to each gene. Default set to \code{"normalized"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of data frames giving CARMA scores per gene.
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
carma.scores <- CARMA(test.data, region.type="arm", hg.version="hg19")
gene.scores <- genesCARMA(carma.scores, "normalized") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
