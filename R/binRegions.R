binRegions <-
function(hg, bin.width, overlap){
	arms		= paste0(hg$chr, "_", hg$arm)
	nreg		= nrow(hg)
	res.reg		= as.data.frame(matrix(NA, nrow=0, ncol=5))
	if(overlap==F){
		for(i in 1:nreg){
			tmp.reg		= hg[i, ]
			larm			= tmp.reg$endpos-tmp.reg$startpos+1
			nbin			= larm/bin.width
			int.bin		= round(nbin)
			int.bin		= ifelse(int.bin<1,1,int.bin)
			lbin			= floor(larm/int.bin)
			bin.tab		= data.frame(chr=tmp.reg[1,1], arm=arms[i], 
							bin=paste0(arms[i], "_", 1:int.bin), 
							startpos=cumsum(c(tmp.reg$startpos, rep(lbin, 
							int.bin-1))), 
							endpos=cumsum(c(tmp.reg$startpos+lbin-1, rep(lbin, 
							int.bin-1))))
			bin.tab$endpos[int.bin] = hg$endpos[i]
			res.reg		= rbind(res.reg, bin.tab)
		}
	} else if(overlap==T){
		for(i in 1:nreg){
			tmp.reg		= hg[i, ]
			larm		= tmp.reg$endpos-tmp.reg$startpos+1
			nbin		= larm/bin.width
			int.bin		= ceiling(nbin)
			lbin		= bin.width-round(((bin.width*int.bin)-larm)/(int.bin-1))
			bin.tab		= data.frame(chr=tmp.reg[1,1], arm=arms[i], 
							bin=paste0(arms[i], "_", 1:int.bin), 
							startpos=cumsum(c(tmp.reg$startpos, 
							rep(lbin, int.bin-1))), 
							endpos=cumsum(c(tmp.reg$startpos+bin.width-1, 
							rep(lbin, int.bin-1))))
			bin.tab$endpos[int.bin] = hg$endpos[i]
			res.reg		= rbind(res.reg, bin.tab)
		}
	}
	return(res.reg)	
}
