carma.S6 <-
function (seg) {
    L = (seg$endpos - seg$startpos)/10^6
	H = (seg$nMajor - seg$nMinor)
	res = sum(L * (H^2))
#    res = sqrt(res)
    	return(res)
}
