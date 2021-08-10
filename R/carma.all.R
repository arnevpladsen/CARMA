carma.all <-
function(seg, ploidy, gender=NULL) {
	res = c()	
	res[length(res)+1] = carma.S1(seg, ploidy, gender) 	# AMP
	res[length(res)+1] = carma.S2(seg, ploidy, gender) 	# DEL
	res[length(res)+1] = carma.S3(seg) 					# STP
	res[length(res)+1] = carma.S4(seg) 					# CRV
	res[length(res)+1] = carma.S5(seg) 					# LOH
	res[length(res)+1] = carma.S6(seg) 					# ASM
	res = round(log2(res+1), 5)
	names(res) = c("AMP", "DEL", "STP", "CRV", "LOH", "ASM")
	return(res)
}
