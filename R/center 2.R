center <-
function(seg) {
    res = c()
    for (i in 1:nrow(seg)) {
    		a = (seg$nMinor + seg$nMajor)[i]
    		b = round((seg$endpos[i]-seg$startpos[i]+1)/10^5)
    		res = c(res, rep(a,b))
    }
    return(round(median(res)))
}
