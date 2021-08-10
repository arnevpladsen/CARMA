#'
#' @export

plotDataset = function(carma, group=NULL, group.order=NULL, plot.titles=NULL, plot.path, n.plot.cols=1, max.color.score=0.7, FUN=median, ...){
	track.gap.x=0.06
	track.gap.y=0.1
	n.groups = ifelse(is.null(group), 1, length(group.order))
	if(is.null(group)){
		group=rep(1, carma[["parameters"]]$number_of_samples)
		group.order = c(1)
	}
	if(is.null(group.order)){
		group.order = unique(group)	
	}
	if(is.null(plot.titles)){
		plot.titles = rep("", n.groups)
	}
	hg.version = carma[["parameters"]][["hg.version"]]
	if(hg.version == "hg19"){
		hg = hg19
	} else if(hg.version == "hg38"){
		hg = hg38
	}
	reg.ref = carma[["start_stop_per_region"]]
	chrom = unique(reg.ref$chr)
	nchr = length(chrom)
	gen.ref = hg[hg$chr %in% chrom, ]
	chr.length = unlist(lapply(chrom, function(x) 
		max(gen.ref$endpos[gen.ref$chr==x])))
	gen.length = sum(chr.length)
	gap = round(gen.length*track.gap.x/nchr)
	indx.gap = gen.length*0.04
	cum.vec = c(0, cumsum((chr.length+gap)[-nchr])) + 
		(indx.gap/2) + gap
	names(cum.vec) = chrom
	gen.ref = cbind(gen.ref, glob.startpos=gen.ref$startpos+
				cum.vec[match(gen.ref$chr, names(cum.vec))], 
				glob.endpos=gen.ref$endpos+ 
				cum.vec[match(gen.ref$chr, names(cum.vec))] )
	reg.ref = cbind(reg.ref, glob.startpos= reg.ref$startpos+ 
		cum.vec[match(reg.ref$chr, names(cum.vec))], 
		glob.endpos=reg.ref$endpos + 
		cum.vec[match(reg.ref$chr, names(cum.vec))])	
	summary.scores = summaryCARMA(carma)

	# Define plotting params for score tracks
	n.plot.rows = ceiling(n.groups/n.plot.cols)
	if(n.plot.rows*n.plot.cols == n.groups){
	n.plot.rows = n.plot.rows+1	
	}
	relative.width = n.plot.cols/n.plot.rows
	height = 4000*n.plot.rows
	col.carma = c(	"#038C7F", "#025951", "#F2CC0F", 
					"#F26E22", "#64BF4B", "#D94929")
	heat.list = heatCol(col.carma)
	png(file= plot.path, width=height*relative.width, height=height)	
	par(mar=rep(0,4))
	plot(1,1, type="n", axes=F, ylab="", xlab="", yaxs="i", 
		xlim=c(-1,2*n.plot.cols-1), ylim=c(-2*n.plot.rows+1,1), 
		xaxs="i", yaxs="i")
	x.offset = rep(seq(0, n.plot.cols*2-1, by=2), n.plot.rows)
	y.offset = rep(seq(0, -n.plot.rows*2, by=-2), each=n.plot.cols)
	for(i in 1: n.groups){
		tmp.group = group.order[i]
		select.samples = which(group == tmp.group) 
		tmp.carma = apply(summary.scores[, select.samples], 1, function(x) 
			FUN(x, na.rm=T))
		carma.scores = vectorToTable(tmp.carma, reg.ref)
		carma.scores = carma.scores[,ncol(carma.scores):1]
		carma.scores[carma.scores>1] = 1
		carma.scores[is.na(carma.scores)] = 0
		score.list = formatPlotScores(carma.scores, reg.ref)
		n.score.tracks = length(score.list)
		y.track.start = 8
		y.max = sum(y.track.start, n.score.tracks, 1)
		x.max = max(gen.ref$glob.endpos) + gap + (indx.gap/2)
		x.step = x.max/1000

		# Score track
		n.scores = length(score.list)
		for(k in 1:n.scores){
			tmp.score = score.list[[k]]
			tmp.score = rbind(tmp.score, tmp.score[1,])
			n.regions = nrow(tmp.score)
			tmp.score$score[n.regions] = 1
			tmp.score$glob.startpos[n.regions] = x.max-
				(indx.gap/2)
			tmp.score$glob.endpos[n.regions] = x.max+
				(indx.gap/2)
			for(m in 1:n.regions){
				reg.score = tmp.score[m,]
				y.start = y.track.start + k - 1
				y.stop = y.start + 1-track.gap.y
				y = c(y.start, y.stop)/y.max
				x0 = reg.score$glob.startpos
				x1 = reg.score$glob.endpos
				x.rect = c(seq(x0, x1, x.step), x1)
				n.x = length(x.rect)
				x.rect = c(x.rect, x.rect[order(x.rect, 
					decreasing=T)])	
				x.circ = sin((x.rect/x.max)*2*pi)*rep(y, 
					each=n.x) + x.offset[i]
				y.circ = cos((x.rect/x.max)*2*pi)*rep(y, 
					each=n.x) + y.offset[i]	
				select.col = round(reg.score$score/max.color.score*100)+1
				select.col = ifelse(select.col > 101, 101, select.col)
				tmp.col = heat.list[[k]][select.col]
				polygon(x=x.circ, y=y.circ, border=NA, 
					col= tmp.col)
				text(0 + x.offset[i],sum(y)/2 + y.offset[i], 
					labels=names(score.list)[k], col="white", font=2, cex=6)
			}
		}
		# Add chromosome numbers
		nchrom = nrow(gen.ref)
		x = ((gen.ref$glob.startpos[seq(1,nchrom,2)] + 
			gen.ref$glob.endpos[seq(2,nchrom,2)])/2)/x.max
		y = (y.track.start-0.3)/y.max
		rot.vec = -360*x
		rot.vec[rot.vec < -90 & rot.vec > -270] = rot.vec[rot.vec < 
			-90 & rot.vec > -270] + 180
		x.circ = sin((x)*2*pi)*y + x.offset[i]
		y.circ = cos((x)*2*pi)*y + y.offset[i]
		n.x = length(x)
		for(j in 1:n.x){
			text(x.circ[j], y.circ[j], labels=chrom[j], col="grey70", 
				font=2, cex=6, srt=rot.vec[j])
		}
		# Add plot title
		text(x.offset[i], y.offset[i], labels=plot.titles[i], cex=11, col="grey50", 
			font=2)
		# Make heatbars
		x.width = 0.7
		x.mid = ((n.groups-((n.plot.rows-1)*n.plot.cols))+n.plot.cols)/2*2-1
		x.start = x.mid-(x.width/2)
		y.mid 	= -(n.plot.rows*2)+2
		y.height = 0.7
		y.gap = 0.02
		indx = colnames(carma[[1]][[1]])[6:1]
		n.indx = length(indx)
		for(j in 1:n.indx){
			tmp.col = heat.list[[j]]
			ncol = length(tmp.col)
			for(k in 1:ncol){
				x.line.start = x.start + (k-1)/ncol*x.width
				x.line.end = x.start + k/ncol*x.width
				y.line.start = y.mid-(y.height/2)+(j-1)*(y.height/n.indx)+(y.gap/2)
				y.line.end = y.line.start + (y.height/n.indx) - y.gap
				rect(x.line.start, y.line.start, x.line.end, 
					y.line.end, col=tmp.col[k], border=NA)
			}
			x.lab = x.line.end + 0.01				
			y.lab = (y.line.start + y.line.end)/2-0.007
			text(x.lab, y.lab, labels=indx[j], 
				col="grey40", cex=7, pos=4, font=2)
		}
		y.grid = y.line.end + 0.005
		segments(x.start, y.grid, x.start+ x.width, lwd=12, 
			lend=1)
		x.grid = x.start + x.width*c(0.0025, 0.5, 0.9975)
		segments(x.grid, y.grid, x.grid, y.grid+0.01, 
			lwd=8, lend=1)
		text(x.grid, rep(y.grid)+0.02, labels=c(0, round(max.color.score/2, 2), 
			max.color.score), cex=6, pos=3)
	}
	dev.off()
}
