\name{carma.S2}
\alias{carma.S2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate CARMA score DEL
}
\description{
This function will not be called directly by the user. The function is used by the \code{carma.all} function to calculate the CARMA score DEL, a measure of the extent of regional loss in copy number.
}
\usage{
carma.S2(seg, ploidy, gender = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a data frame with copy number data corresponing to one region.
}
  \item{ploidy}{
an integer value corresponding to the median copy number across the whole genome, calculated using the \code{center} function.
}
  \item{gender}{
a character string giving gender. If \code{gender="male"}, the function will correct for an expected lower value for the X chromosome in segments covering the X chromosome. This is done by dividing the ploidy by 2 when calculating AMP and DEL in these segments.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric score.
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
seg <- segment[segment$chr==1, ] 
ploidy <- center(test.data[[1]])
carma.S2(seg, ploidy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
