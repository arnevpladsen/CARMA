vectorToTable = function(tmp.carma, reg.ref){
	indx = unique(substr(names(tmp.carma), 1,3))
	n.indx = length(indx)
	regs = reg.ref$region
	n.regs = length(regs)
	res.mat = as.data.frame(matrix(NA, n.regs, n.indx, dimnames=list(regs, indx)))
	tmp.carma[tmp.carma==Inf] = NA
	for(i in 1:n.indx){
		for(j in 1:n.regs){
			select.score = which( substr(names(tmp.carma), 1,3) == indx[i] &
				 substr(names(tmp.carma), 5,nchar(names(tmp.carma))) == regs[j] )
			if(length(select.score) == 0){next}
			res.mat[j,i] = tmp.carma[select.score]
		}	
	}
	res.mat = round(res.mat, 3)
	return(res.mat)	
}