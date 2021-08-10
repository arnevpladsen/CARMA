\name{center}
\alias{center}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate median ploidy across the genome
}
\description{
This function will not be called directly by the user. The function is used by the \code{CARMA} function and takes as argument a segment table for a region and returns the median copy number in that region. 
}
\usage{
center(seg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a data frame with copy number data for a region.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric value
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
segment <- test.data[[1]]
center(segment)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
