#'
#' @export


CARMA <-
function(segment.list, region.type=c("arm", "chrom", "bin") , bin.width=5*10^7, overlap=T, hg.version=c("hg38","hg19"), gender=NULL, exclude.x.chrom=T, normalization.set=NULL){
	# Format input segments
	reg.list = vector("list", length(segment.list))
	names(reg.list) = names(segment.list)
	seg.list = reg.list
	hg.version = hg.version[1]
	region.type = region.type[1]
	if(hg.version=="hg19"){
		hg	= hg19
	} else if(hg.version=="hg38"){
		hg = hg38
	}
	bin.regions = NULL
	if(region.type == "bin"){
		bin.regions	= binRegions(hg, bin.width, overlap)	
	}
	for(i in 1:length(segment.list)){	
		seg	= segment.list[[i]]
		seg	= continuousSegments(seg, hg, exclude.x.chrom)
		seg.list[[i]] = seg
		if(region.type %in% c("arm", "chrom")){
			reg.list[[i]] = seg
		} else if(region.type == "bin"){
			seg = binSegments(seg, bin.regions)
			reg.list[[i]] = seg	
		}
	}
	# Run CARMA algorithm
	carma.res.reg = vector("list", length(reg.list))
	names(carma.res.reg) = names(reg.list)
	carma.res.genome	 = carma.res.reg
	for(i in 1:length(reg.list)){
		seg	= reg.list[[i]]
		ploidy = center(seg)
		indx	 = names(carma.all(seg[1,], ploidy, NULL))
		if(region.type %in% c("arm", "chrom")){
			arms = unique(seg$arm) # calculate armwise first
			res.reg	= as.data.frame(matrix(NA, nrow=length(arms), ncol=length(indx), 
						dimnames=list(arms, indx)))	
			for(j in 1:length(arms)){
				tmp.seg	= seg[which(seg$arm == arms[j]), ]
				res.reg[j,]	= carma.all(tmp.seg, ploidy, gender[i])
			}
			if(region.type == "chrom"){
				res.reg	= reformatCARMA(res.reg, "chrom")
			} 
			carma.res.reg[[i]] = res.reg
		} else if(region.type=="bin"){
			bins = unique(seg$bin)
			res.reg	= as.data.frame(matrix(NA, nrow=length(bins), ncol=length(indx), 
						dimnames=list(bins, indx)))
			for(j in 1:length(bins)){
				tmp.seg = seg[which(seg$bin == bins[j]), ]
				res.reg[j,]	= carma.all(tmp.seg, ploidy, gender[i])
			}
			carma.res.reg[[i]] = res.reg
		}
		carma.res.genome[[i]] = reformatCARMA(res.reg, "whole_genome")
	}
	# Normalize CARMA scores
	carma.norm.res = normalizeCARMA(carma.res.reg, carma.res.genome, 
		indx, region.type, bin.regions, exclude.x.chrom, normalization.set, 
		quant=0.99)
	# Return result file
	region.positions = getReference(region.type, hg, bin.width, overlap, 
		exclude.x.chrom)
	params = list(number_of_samples=length(segment.list), region.type=region.type, 
		bin.width=bin.width, overlap=overlap, hg.version=hg.version, gender=gender, 
		exclude.x.chrom=exclude.x.chrom, 
		normalization.set=ifelse(is.null(normalization.set), 
		"NULL", "external normalization set used"))
	res=list(	raw_regional_CARMA_scores=carma.res.reg, 
				raw_genome_wide_CARMA_scores=carma.res.genome,
				normalized_regional_CARMA_scores= carma.norm.res[["regional"]],
				normalized_genome_wide_CARMA_scores= carma.norm.res[["genome"]],
				start_stop_per_region = region.positions,
				parameters = params,
				segments = seg.list
				)
	return(res)
}
