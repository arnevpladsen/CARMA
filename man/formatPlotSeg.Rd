\name{formatPlotSeg}
\alias{formatPlotSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format copy number segments for use in plotting functions
}
\description{
This function will not be called directly by the user. The function is used by the \code{plotSample} function to reformat copy number segments to an appropriate format.
}
\usage{
formatPlotSeg(seg, cum.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a table with copy number segments obtained from the output file produced by the \code{CARMA} algorithm, found in the \code{[["segments"]]} list.
}
  \item{cum.vec}{
a vector with global positions for each chromosome.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame giving global start and end positions for each copy number segment.
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
