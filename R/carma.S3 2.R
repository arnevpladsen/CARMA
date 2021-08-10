carma.S3 <-
function(seg) {
	seg.start = seg$startpos
	seg.stop = seg$endpos
	seg.mean = seg$nMinor + seg$nMajor
	# 1st derivative
	x = (seg.start+seg.stop)/2
	y = seg.mean
	seg.mean = 0
    if (length(x)>1) {
    	k = 1:(length(x)-1)
        dx = 1/tanh(1/((x[k+1]-x[k])/10^6))
        dy = tanh(y[k+1]-y[k])
        seg.start = x[k]
        seg.stop = x[k+1]
        seg.mean = dy/dx
    }
    # area under curve
    L = (seg.stop-seg.start)/10^6
    res = sum(L * seg.mean^2)
#    res = sqrt(res)
    return(res)	
}
