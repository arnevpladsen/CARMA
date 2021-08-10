formatPlotSeg = function(seg, cum.vec){
	glob.start = seg$startpos+ cum.vec[match(seg$chr, names(cum.vec))]
	glob.end = seg$endpos+ cum.vec[match(seg$chr, names(cum.vec))]	
	mean = seg$nMajor + seg$nMinor
	res = data.frame(chr=seg$chr, glob.startpos=glob.start, glob.endpos=glob.end, 
		mean=mean)
	return(res)
}