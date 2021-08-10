carma.S5 <-
function(seg) {
	L = (seg$endpos - seg$startpos)/10^6
	res = sum(L[seg$nMinor == 0]) 
    	return(res)
}
