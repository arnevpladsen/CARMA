\name{combineTables}
\alias{combineTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine CARMA score tables
}
\description{
This function will not be called directly by the user. The function is used by the \code{normalizeCARMA} function to combine multiple CARMA score tables into one.
}
\usage{
combineTables(regs, carma.res.reg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regs}{
a vector giving the region ID.
}
  \item{carma.res.reg}{
a list of length corresponding to the number of samples, each containing regional carma scores for a given sample.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing regional CARMA scores for all samples. Rows represent regions, while columns represent samples.
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
