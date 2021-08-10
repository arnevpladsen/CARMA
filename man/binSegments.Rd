\name{binSegments}
\alias{binSegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split a segmented copy number profile into a set of predefined bins
}
\description{
This function will not be called directly by the user. The function splits a segmented copy number profile into a set of predefined bins, defined by the \code{binRegions} function.
}
\usage{
binSegments(seg, bin.regions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a data frame containing the input copy number profile. Output from the \code{continuousSegments} function.
}
  \item{bin.regions}{
a data frame with predefinded bin regions. Output from the \code{binRegions} function.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
seg = continuousSegments(seg=segment,hg=hg19, exclude.x.chrom=F)
bin.reg <- binRegions(hg=hg19, bin.width=3*10^7, overlap=T)
bin.seg <- binSegments(seg=seg, bin.regions=bin.reg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
