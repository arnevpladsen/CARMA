getReference <-
function(region.type, hg, bin.width, overlap, exclude.x.chrom){
	if(region.type=="arm"){
		ref.region = data.frame(region=paste0(hg$chr, "_", hg$arm), chr=hg$chr,
			startpos=hg$startpos, endpos=hg$endpos)
	} else if(region.type=="chrom"){
		ref.region = data.frame(region=hg$chr[seq(1,46,2)],
								chr=hg$chr[seq(1,46,2)],
								startpos=hg$startpos[seq(1,46,2)],
								endpos=hg$endpos[seq(2,46,2)])
	} else if(region.type=="bin"){
		bin.reg = binRegions(hg, bin.width, overlap)
		ref.region = data.frame(region=bin.reg$bin, chr=bin.reg$chr, 
			startpos=bin.reg$startpos, endpos=bin.reg$endpos)	
	}
	if(exclude.x.chrom == TRUE){
		ref.region = ref.region[which(!ref.region$chr == "X"), ]
	}
	return(ref.region)
}
