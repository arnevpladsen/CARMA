splitTable <-
function(norm.res, carma.res.reg){
	carma.norm.res.reg = vector("list", length(carma.res.reg))
	names(carma.norm.res.reg) = names(carma.res.reg)
	indx = colnames(carma.res.reg[[1]])
	for(i in 1:length(carma.norm.res.reg)){
		tmp.vec = norm.res[,i]
		tmp.tab = matrix(tmp.vec, ncol=length(indx))
		colnames(tmp.tab) = indx
		reg.names = names(tmp.vec)[1:nrow(tmp.tab)]
		rownames(tmp.tab) = substr(reg.names, 5, nchar(reg.names))
		exclude = which(apply(tmp.tab, 1, function(x) length(which(is.na(x)))) ==
					length(indx))
		if(length(exclude)>0){
			tmp.tab = tmp.tab[-exclude, ]
		}
		carma.norm.res.reg[[i]] = tmp.tab
	}	
	return(carma.norm.res.reg)
}
