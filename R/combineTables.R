combineTables <-
function(regs, carma.res.reg){
	res.tab = as.data.frame(matrix(NA, nrow=length(regs), 
		ncol=length(carma.res.reg), dimnames=list(regs, names(carma.res.reg)) )) 
	nsamp = length(carma.res.reg)
	for(i in 1:nsamp){
		tmp.res = carma.res.reg[[i]]
		tmp.vec = unlist(tmp.res)
		names(tmp.vec) = paste0(rep(colnames(tmp.res), each=nrow(tmp.res)), "_",
							rownames(tmp.res))
		res.tab[match(names(tmp.vec), rownames(res.tab)), i] = tmp.vec
	}
	return(res.tab)
}
