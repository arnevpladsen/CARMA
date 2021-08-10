# CARMA

## Description

The package provides functions for calculating CARMA scores within genomic regions from segmented DNA copy number profiles. The CARMA scores are numeric scores that reflect the presence of specific genomic motifs within the selected regions.

## Installation
You can install the current CARMA version using the command
```R
devtools::install_github("arnevpladsen/CARMA")
```

## Details

The main function in the carma package is the ```CARMA``` function. Input for the function is a list of length samples, containing data frames with segmented integer copy number data. Test data is supplied in the object ```test.data``` and shows the format of the input data. CARMA scores can be converted to summary tables using the ```summaryCARMA``` function. Regional CARMA scores can also be mapped to genes using the ```genesCARMA``` function. The package addtionally supplies plotting functions: per sample using the ```plotSample``` function, or across a dataset using the ```plotDataset``` function.

## Authors

Authors: Arne Pladsen [aut, cre], Gro Nilsen [aut], Ole Christian Lingjærde [aut]

Maintainer: Arne Pladsen <arnpla@rr-research.no>

## References

https://www.nature.com/articles/s42003-020-0884-6

## Examples
```R
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
```