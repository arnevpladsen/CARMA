continuousSegments <-
function(seg, hg, exclude.x.chrom=F) {
	seg = seg[which(!seg$chr == "Y"), ]
	if(exclude.x.chrom==T){
		seg = seg[which(!seg$chr == "X"), ]	
	}
	narms = nrow(hg)
  	res = seg[0, ]
  	for(i in seq(1,narms, 2)){
   		tmp.seg		= seg[which(seg$chr == hg$chr[i]), ]
   		nseg			= nrow(tmp.seg)
   		if(nseg > 0){
			# Check for segment spanning the centromere
    		select.span	= which(tmp.seg$startpos < hg$endpos[i] & 
    						tmp.seg$endpos > hg$startpos[i+1])
	    		if(length(select.span) == 1){
	    			span.seg	= rbind(tmp.seg[select.span,], 
	    				tmp.seg[select.span,])
	    			span.seg$endpos[1]		= hg$endpos[i]
	    			span.seg$startpos[2]	= hg$startpos[i+1]
				split.seg	= tmp.seg[0,]		
	    			if(select.span > 1){	
	    				split.seg	= rbind(tmp.seg[1:(select.span-1), ], 
	    					span.seg)
	    			} else {
	    				split.seg	= span.seg
	    			}
	    			if(select.span < nseg){
	    				split.seg	= rbind(split.seg, tmp.seg[
	    					(select.span +1): nseg, ])
	    			} 
	    			tmp.seg	= split.seg	
	    		}
	    		# Fill gap between segments for each chromosome arm
	    		gap.seg	= tmp.seg[0,]
	    		for(j in 0:1){
	    			arm.seg	= tmp.seg[which(tmp.seg$startpos < hg$endpos[i+j] &
	    						tmp.seg$startpos >= hg$startpos[i+j]), ]
	    			
	    			if(nrow(arm.seg) == 0){
	    				next
	    			} else if(nrow(arm.seg) == 1){
	    				arm.seg$arm	= paste0(hg$chr[i], "_", c("p","q")[j+1])
	    				gap.seg 		= rbind(gap.seg, arm.seg)
	    				
	    			} else {
	    				arm.seg$arm	= paste0(hg$chr[i], "_", c("p","q")[j+1])
	    				nseg 	= nrow(arm.seg)
	    				gap		= arm.seg$startpos[2:nseg] - 
	    							arm.seg$endpos[1:(nseg-1)] - 1
	    				half	= round(gap/2)
	    				arm.seg$endpos[1:(nseg-1)] = arm.seg$endpos[1:(nseg-1)]+ 
	    												half
	    				arm.seg$startpos[2:nseg] = arm.seg$endpos[1:(nseg-1)] + 1 
	    				gap.seg = rbind(gap.seg, arm.seg)		
	    			}
	    			
	    		}    		   		
			res	= rbind(res, gap.seg)
		}		
	}
	res	= data.frame(chr=res$chr, arm=as.character(res$arm), startpos=res$startpos, 
		endpos=res$endpos, nMajor=res$nMajor, nMinor=res$nMinor)
	return(res)	
}
