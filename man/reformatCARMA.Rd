\name{reformatCARMA}
\alias{reformatCARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reformat CARMA scores
}
\description{
This function will not be called directly by the user. The function is used by the \code{CARMA} function and converts a CARMA results matrix with armwise CARMA scores per row into a matrix where each row represents scores per whole chromosome or whole genome.
}
\usage{
reformatCARMA(res.reg, region.type = "chrom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.reg}{
a CARMA results matrix
}
  \item{region.type}{
the defined region type. Either \code{"chrom"} or \code{"whole_genome"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A CARMA results matrix where each row represents scores per whole chromosome or whole genome
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
