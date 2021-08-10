\name{carma-package}
\alias{carma-package}
\alias{carma}
\docType{package}
\title{
\packageTitle{carma}
}
\description{
\packageDescription{carma}
}


\details{
%~~ The DESCRIPTION file:
%\packageDESCRIPTION{carma}
%\packageIndices{carma}
The main function in the \code{carma} package is the \code{CARMA} function. Input for the function is a list of length samples, containing data frames with segmented integer copy number data. Test data is supplied in the object \code{test.data} and shows the format of the input data. CARMA scores can be converted to summary tables using the \code{summaryCARMA} function. Regional CARMA scores can also be mapped to genes using the \code{genesCARMA} function. The package addtionally supplies plotting functions: per sample using the \code{plotSample} function, or across a dataset using the \code{plotDataset} function. 
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}

\author{
Authors: \packageAuthor{carma}

Maintainer: \packageMaintainer{carma}
}

\references{
http://biorxiv.org/lookup/doi/10.1101/769356
%~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# Format of input copy number profiles
head(test.data[[1]])

# Calculate CARMA scores per chromosome arm
carma.res <- CARMA(test.data, region.type="arm", hg="hg19")

# Calculate CARMA scores per whole chromosome
carma.res <- CARMA(test.data, region.type="chrom", hg="hg19")

# Calculate CARMA scores within 3*10^7 base pairs long regions
carma.res <- CARMA(test.data, region.type="bin", bin.width=3*10^7, hg="hg19")

# Plot CARMA scores per sample 
carma.res <- CARMA(test.data, region.type="arm", hg="hg19")
plotSample(carma.res, samples=c(1,2,3), plot.dir="~/Desktop/")

# Plot CARMA scores across the whole dataset
plotDataset(carma.res, plot.path="~/Desktop/carma_scores_entire_dataset.png")

# Plot CARMA scores across groups in the dataset
plotDataset(carma.res, group=rep(c("pos", "neg"), 50), 
	group.order=c("neg", "pos"), plot.titles=c("Negative", "Positive"), 
	plot.path="~/Desktop/carma_scores_groups.png", n.plot.cols=3)
}