carma.S2 <-
function(seg, ploidy, gender=NULL) {
  	if(!is.null(gender)){
 		if(gender=="male" & seg$chr[1]=="X"){
 			ploidy	= ploidy/2
 		}
 	}
    L = (seg$endpos - seg$startpos)/10^6
    H = (seg$nMinor + seg$nMajor) - ploidy
    res = sum(L[H < 0] * (H[H < 0])^2)
#    res = sqrt(res)
    return(res)
}
