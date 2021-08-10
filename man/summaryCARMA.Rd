\name{summaryCARMA}
\alias{summaryCARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make summary table of CARMA scores for all samples
}
\description{
This function reformats the output data from the \code{CARMA} function into one single table containing CARMA scores for all samples. 
}
\usage{
summaryCARMA(carma, scores = "regions", data.type = "normalized")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{carma}{
the output file from the \code{CARMA} function, containing CARMA scores.
}
  \item{scores}{
the type of scores to be returned: \code{"regions"} or \code{"genome_wide"}, corresponding to regional or aggregated genome wide scores. Default set to \code{"regions"}.
}
  \item{data.type}{
the type of data to be returned: \code{"normalized"} or \code{"raw"}. Default set to \code{"normalized"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing CARMA scores for all samples.
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
summary.table <- summaryCARMA(carma.scores) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
