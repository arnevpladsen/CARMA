#'
#' @export
#'

genesCARMA <-
function(carma, data.type=c("normalized", "raw")){
	data.type = data.type[1]
	if(data.type=="normalized"){
		score.list = carma[["normalized_regional_CARMA_scores"]]		
	} else if(data.type=="raw"){
		score.list = carma[["raw_regional_CARMA_scores"]]				
	}
	if(carma[["parameters"]][["hg.version"]] == "hg19"){
		genes.tab = hg19_genes
	} else if(carma[["parameters"]][["hg.version"]] == "hg38"){
		genes.tab = hg38_genes
	}
	ref.pos = carma[["start_stop_per_region"]]
	nsamp = length(score.list)
	res.list = vector("list", nsamp)
	names(res.list) = names(score.list)
	for(i in 1:nsamp){
		tmp.res = score.list[[i]]
		nreg = nrow(tmp.res)
		indx = colnames(tmp.res)
		tmp.gene = cbind(genes.tab, as.data.frame(matrix(NA, nrow=nrow(genes.tab), 
					ncol=length(indx), dimnames=list(1:nrow(genes.tab), indx))) )
		for(j in 1:nreg){
					select.ref = match(rownames(tmp.res)[j], ref.pos$region)
					tmp.chr = ref.pos$chr[select.ref]
					tmp.start = ref.pos$startpos[select.ref] 
					tmp.end = ref.pos$endpos[select.ref] 
					select.row = which(	tmp.gene$chr %in% tmp.chr &
									 	tmp.gene$startpos >= tmp.start &
									 	tmp.gene$endpos <= tmp.end)
					select.col = which(colnames(tmp.gene) %in% indx)
					tmp.gene[select.row, select.col] = rep(tmp.res[j, ], 
						each=length(select.row))						
		}
		res.list[[i]] = tmp.gene
	}
	return(res.list)	
}
