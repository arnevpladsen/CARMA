\name{vectorToTable}
\alias{vectorToTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format CARMA scores for use in plotting function
}
\description{
This function will not be called directly by the user. The function is used by the \code{plotDataset} function to reformat CARMA scores to an appropriate format.
}
\usage{
vectorToTable(tmp.carma, reg.ref)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmp.carma}{
a vector giving region specific CARMA scores.
}
  \item{reg.ref}{
a table giving chromosome as well as start and stop positions for each CARMA region.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing CARMA scores.
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
