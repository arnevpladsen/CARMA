\name{continuousSegments}
\alias{continuousSegments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make continuous segments and split segments spanning the centromere
}
\description{
This function will not be called directly by the user. The function is used by the \code{CARMA} function to reformat input copy number profiles by extending segments flanked by gaps and by splitting segments spanning the centromere. The function can also optionally exclude the X chromosome.
}
\usage{
continuousSegments(seg, hg, exclude.x.chrom = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
a data frame containing the copy number profile of a single sample.
}
  \item{hg}{
a data frame containing start and end positions of chromosome arms: \code{hg19} or \code{hg38}.
}
  \item{exclude.x.chrom}{
if \code{exclude.x.chrom=TRUE}, the X chromosome is excluded in the output file.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing formatted segment start and end positions.
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
continuousSegments(segment, hg19, exclude.x.chrom=F)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
