#'
#' @export
#'

summaryCARMA <-
function(carma, scores=c("regions", "genome_wide"), data.type=c("normalized", "raw")){
	scores = scores[1]
	data.type = data.type[1]
	if(scores == "regions"){
		if(data.type=="normalized"){
			carma.res = carma[["normalized_regional_CARMA_scores"]]
		} else if(data.type=="raw"){
			carma.res = carma[["raw_regional_CARMA_scores"]]
		}
	} else if(scores == "genome_wide"){
		if(data.type=="normalized"){
			carma.res = carma[["normalized_genome_wide_CARMA_scores"]]
			carma.res = lapply(carma.res, function(x) as.data.frame(matrix(x, 
				nrow=1, dimnames=list("", names(x)) )))
		} else if(data.type=="raw"){
			carma.res = carma[["raw_genome_wide_CARMA_scores"]]
			carma.res = lapply(carma.res, function(x) as.data.frame(matrix(x, 
				nrow=1, dimnames=list("", names(x)) )))
		}
	}
	nsamp = length(carma.res)
	indx = colnames(carma.res[[1]])
	ref.tab = carma[["start_stop_per_region"]]
	regs = ref.tab$region
	if(scores=="regions"){
		reg.names = paste0(rep(indx, each=length(regs)), "_", regs)
	} else if(scores=="genome_wide"){
		reg.names = indx
	}
	summary.tab = as.data.frame(matrix(NA, nrow=length(reg.names), 
		ncol=nsamp, dimnames=list(reg.names, names(carma.res)) ))
	for(i in 1:nsamp){
		tmp.tab = as.data.frame(carma.res[[i]])
		tmp.vec = unlist(tmp.tab)
		if(scores=="regions"){
		names(tmp.vec) = paste0(rep(colnames(tmp.tab), each=nrow(tmp.tab)), 
			"_", rownames(tmp.tab))
		} else if(scores=="genome_wide"){
			names(tmp.vec) = indx
		}
		summary.tab[match(names(tmp.vec), rownames(summary.tab)),i] = tmp.vec		
	}
	exclude = which(apply(summary.tab, 1, function(x) length(which(is.na(x)))) == 
		ncol(summary.tab)  )
	if(length(exclude) > 0){
		summary.tab = summary.tab[-exclude, ]
	}
	return(summary.tab)
}
