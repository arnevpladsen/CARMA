formatPlotScores = function(scores, reg.ref){
	n.scores = ncol(scores)
	score.names = colnames(scores)
	res.list = vector("list", n.scores)
	names(res.list) = score.names
	for(i in 1:n.scores){
		tmp.score = cbind(reg.ref, score=0)
		tmp.score$score[match(rownames(scores),  tmp.score$region)] = scores[,i]
		res.list[[i]] = tmp.score
	}
	return(res.list)
}