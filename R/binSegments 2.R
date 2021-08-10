binSegments <-
function(seg, bin.regions) {
	nbins		= nrow(bin.regions)
	res.seg		= data.frame(matrix(NA, nrow=0, ncol=7))
	for(i in 1:nbins){
		tmp.bin		= bin.regions[i, ]
		match.seg	= which(as.character(seg$chr)==as.character(tmp.bin$chr) &
						seg$startpos<tmp.bin$endpos & seg$endpos>tmp.bin$startpos)
		if(length(match.seg) == 0){next}
		tmp.seg		= seg[match.seg, ]
		nseg		= nrow(tmp.seg)
		tmp.seg$startpos[1] = ifelse(tmp.seg$startpos[1] < tmp.bin$startpos, 
								tmp.bin$startpos, tmp.seg$startpos[1])
		tmp.seg$endpos[nseg] = ifelse(tmp.seg$endpos[nseg]>tmp.bin$endpos, 
								tmp.bin$endpos, tmp.seg$endpos[nseg])
		
		tmp.seg		= data.frame(chr=tmp.seg$chr, arm=tmp.seg$arm, bin=tmp.bin$bin,
						startpos=tmp.seg$startpos, endpos=tmp.seg$endpos, 
						nMajor=tmp.seg$nMajor, nMinor=tmp.seg$nMinor)
		res.seg		= rbind(res.seg, tmp.seg)
	}
	return(res.seg)
}
