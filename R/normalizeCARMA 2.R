normalizeCARMA <-
function(carma.res.reg, carma.res.genome, indx, region.type="arm", 
	bin.regions=NULL, exclude.x.chrom=T, normalization.set=NULL, quant=0.99){
		# Normalize regional scores
		if(region.type=="arm"){
			regs = paste0( rep(c(1:22,"X"), each=2), "_", rep(c("p","q"), 23))
		} else if(region.type=="chrom"){
			regs = c(1:22,"X")	
		} else if(region.type=="bin"){
			regs = bin.regions$bin		
		}
		if(exclude.x.chrom==T){
			regs = regs[which(!substr(regs,1,1) == "X")]
		}
		regs = paste0(rep(indx, each=length(regs)), "_", regs)
		comb.res = combineTables(regs, carma.res.reg)
		if(is.null(normalization.set)){
			norm.vector = apply(comb.res, 1, function(x) quantile(x, quant, 
				na.rm=T))
		} else {
			norm.set = normalization.set[["raw_regional_CARMA_scores"]]
			comb.norm = combineTables(regs, norm.set)
			norm.vector = apply(comb.norm, 1, function(x) quantile(x, quant, 
				na.rm=T))	
		}
		norm.res = apply(comb.res, 2, function(x) round(x/norm.vector, 3))
		carma.norm.res.reg = splitTable(norm.res, carma.res.reg)
		# Normalize whole genome scores
		if(is.null(normalization.set)){
			norm.tab = matrix(unlist(carma.res.genome), ncol=length(indx), 
				dimnames=list(names(carma.res.genome), indx), byrow=T)
			norm.vector = apply(norm.tab, 2, function(x) quantile(x, quant, 
				na.rm=T))
		} else {	
			norm.set = normalization.set[["raw_genome_wide_CARMA_scores"]]
			norm.tab = matrix(unlist(norm.set), ncol=length(indx), 
				dimnames=list(names(norm.set), indx), byrow=T)
			norm.vector = apply(norm.tab, 2, function(x) quantile(x, quant, 
				na.rm=T))
		}
		carma.norm.res.genome = lapply(carma.res.genome, function(x) round( 
			x/norm.vector, 3))
		carma.norm.res = list(	regional=carma.norm.res.reg, 
								genome= carma.norm.res.genome)
		return(carma.norm.res)	
}
