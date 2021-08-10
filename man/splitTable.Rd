\name{splitTable}
\alias{splitTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split a CARMA score table
}
\description{
This function will not be called directly by the user. The function is used by the \code{normalizeCARMA} function to split a CARMA score table into a list containing sample-wise CARMA scores.

}
\usage{
splitTable(norm.res, carma.res.reg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{norm.res}{
a data frame containing normalized regional carma scores for all samples.
}
  \item{carma.res.reg}{
a list of length corresponding to the number of samples, each containing raw regional carma scores for a given sample.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of length corresponding to the number of samples, each containing normalized regional carma scores for a given sample.
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
