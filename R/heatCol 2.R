heatCol = function(col.vec){
	n.col = length(col.vec)
	res.list = list()	
	for(i in 1:n.col){
		tmp.col = col.vec[i]	
		rgb.col = col2rgb(tmp.col)
		r = seq(247/255, rgb.col[1]/255, length = 1000)
		b = seq(247/255, rgb.col[2]/255, length = 1000)
		g = seq(247/255, rgb.col[3]/255, length = 1000)			
		x = round((seq(0, 1, length=101)^3)*1000)
		x[x==0] = 1
		res.list[[i]] = rgb(r[x], b[x], g[x])			
	}
	return(res.list)	
}