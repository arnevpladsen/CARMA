\name{carma.S6}
\alias{carma.S6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate CARMA score ASM
}
\description{
This function will not be called directly by the user. The function is used by the \code{carma.all} function to calculate the CARMA score ASM, a measure of the extent of regional assymetry.
}
\usage{
carma.S6(seg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a data frame with copy number data corresponing to one region.
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
carma.S6(seg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
