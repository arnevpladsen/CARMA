reformatCARMA <-
function(res.reg, region.type="chrom"){
	if(region.type=="chrom"){	
		chr	= substr(rownames(res.reg), 1, unlist(gregexpr("_", rownames(res.reg)))-1)
		res = res.reg[0,]
		for(i in c(1:22,"X","Y")){
			tmp.res	= res.reg[which(chr == i), ]
			if(nrow(tmp.res) == 0){next}
			res = rbind(res, apply(tmp.res, 2, function(x) log2(sum(2^x-1, na.rm=T)+1)))
			rownames(res)[nrow(res)] = i
		}
		colnames(res) = colnames(res.reg)
		return(res)
	} else if(region.type=="whole_genome"){
		res = apply(res.reg, 2, function(x) log2(sum(2^x-1, na.rm=T)+1))
		return(res)	
	}
}
