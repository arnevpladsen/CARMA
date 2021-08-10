#'
#' @export

plotSample = function(carma, samples=c(1), plot.dir="~/Desktop/"){
	track.gap.x=0.06
	track.gap.y=0.1
	module.gap=0.2
	nsamples = length(samples)
	hg.version = carma[["parameters"]][["hg.version"]]
	if(hg.version == "hg19"){
		hg = hg19
	} else if(hg.version == "hg38"){
		hg = hg38
	}
	reg.ref = carma[["start_stop_per_region"]]
	for(i in samples){
		samp.name = ifelse(class(i)=="character", i, 
			names(carma[[1]][i]))
		relative.width = 1.3
		height = 4000
		png(file=paste0(plot.dir,"Circos_plot_", samp.name, ".png"), 
			width= height* relative.width, height= height)
			scores = carma[["normalized_regional_CARMA_scores"]][[i]]
			scores = scores[,ncol(scores):1]
			seg =  carma[["segments"]][[i]]
			scores[scores>1] = 1
			scores[is.na(scores)] = 0
			chrom = unique(formatPlotScores(scores, reg.ref)[[1]]$chr)
			nchr = length(chrom)
			gen.ref = hg[hg$chr %in% chrom, ]
			chr.length = unlist(lapply(chrom, function(x) 
				max(gen.ref$endpos[gen.ref$chr==x])))
			gen.length = sum(chr.length)
			gap = round(gen.length*track.gap.x/nchr)
			indx.gap = gen.length*0.03
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
			score.list = list()
			score.list[[1]] = formatPlotScores(scores, reg.ref)
			n.modules = length(score.list)
			n.score.tracks = lapply(score.list, function(x) length(x))
			# Define plotting params for score tracks
			y.track.start = 8
			y.max = sum(y.track.start, unlist(n.score.tracks), 
				module.gap*n.modules, 3)
			x.max = max(gen.ref$glob.endpos) + gap + (indx.gap/2)
			x.step = x.max/1000
			carma.col = c(	"#038C7F", "#025951", "#F2CC0F", 
							"#F26E22", "#64BF4B", "#D94929")
			heat.list = heatCol(carma.col)
			indx.names = list()
			indx.names[[1]] = colnames(scores)
	
			# Define plotting params for copy number track
			ploidy = log2(center(seg)+1)
			seg = formatPlotSeg(seg, cum.vec)
			seg$mean[seg$mean>15] = 16
			cn.max = ifelse(max(seg$mean, na.rm=T) == 16, 16, 
				max(seg$mean, na.rm=T)+1)
			seg$mean = log2(seg$mean +1)
			cn.max = log2(cn.max+1)
			seg$y = seg$mean/cn.max
			ploidy.y = ploidy/cn.max 
			height.cn.track = 2
			
			# Start plotting		
			par(mar=rep(0,4))
			plot(1,1, type="n", axes=F, ylab="", xlab="", yaxs="i", 
				xlim=c(-1,2*relative.width-1), ylim=c(-1,1), 
				xaxs="i", yaxs="i")
				# Score track
				for(j in 1:n.modules){
					tmp.score.list = score.list[[j]]
					n.scores = length(tmp.score.list)
					for(k in 1:n.scores){
						tmp.score = tmp.score.list[[k]]
						tmp.score = rbind(tmp.score, tmp.score[1,])
						n.regions = nrow(tmp.score)
						tmp.score$score[n.regions] = 1
						tmp.score$glob.startpos[n.regions] = x.max-
							(indx.gap/2)
						tmp.score$glob.endpos[n.regions] = x.max+(indx.gap/2)
						for(m in 1:n.regions){
							reg.score = tmp.score[m,]
							y.start = y.track.start + (j-1)*module.gap + (k-1)
							y.stop = y.start + 1-track.gap.y
							y = c(y.start, y.stop)/y.max
							x0 = reg.score$glob.startpos
							x1 = reg.score$glob.endpos
							x.rect = c(seq(x0, x1, x.step), x1)
							n.x = length(x.rect)
							x.rect = c(x.rect, x.rect[order(x.rect, 
								decreasing=T)])	
							x.circ = sin((x.rect/x.max)*2*pi)*rep(y, each=n.x)
							y.circ = cos((x.rect/x.max)*2*pi)*rep(y, each=n.x)	
							select.col = round(reg.score$score*100)+1
							tmp.col = heat.list[[k]][select.col]
							polygon(x=x.circ, y=y.circ, border=NA, 
								col= tmp.col)
							text(0,sum(y)/2, labels=indx.names[[j]][k], 
								col="white", font=2, cex=5)
						}
					}
				}
				# Copy number track
				for(j in 1:nchr){
					x0 = min(gen.ref$glob.startpos[gen.ref$chr %in% chrom[j]])
					x1 = max(gen.ref$glob.endpos[gen.ref$chr %in% chrom[j]])
					x = c(seq(x0, x1, x.step), x1)
					n.x = length(x)
					x = c(x, x[order(x, decreasing=T)])
					y.start = y.track.start + sum(unlist(n.score.tracks)) + 
						module.gap*n.modules
					y.end = y.start+height.cn.track
					y = c(y.start, y.end)/y.max
					y = rep(y, each=n.x)				
					x.circ = sin((x/x.max)*2*pi)*y
					y.circ = cos((x/x.max)*2*pi)*y
					polygon(x=x.circ, y=y.circ, border="grey50", 
						col= "grey95", lwd=0.5)
					y.pl = (y.start + ploidy.y*height.cn.track)/y.max
					x = x[1:n.x]
					x.pl = sin((x/x.max)*2*pi)*y.pl
					y.pl = cos((x/x.max)*2*pi)*y.pl
					lines(x.pl, y.pl, col="grey80", lwd=7, lend=1)
					s = seg[seg$chr %in% chrom[j], ]
					x0.seg = s$glob.startpos
					x1.seg = s$glob.endpos
					x.seg = data.frame(x0.seg, x1.seg)
					n.seg = nrow(x.seg)
					if(n.seg>1){
					x.seg[,2][1:(n.seg-1)] = floor((x.seg[,2][1:(n.seg-1)]
						+ x.seg[,1][2:n.seg])/2)
					x.seg[,1][2:n.seg] = x.seg[,2][1:(n.seg-1)]+1
					
					x.seg.list = apply(x.seg, 1, function(x) c(seq(x[1], 
						x[2], x.step), x[2]))
					} else {
						x.seg.list = vector("list", 1)
						x.seg.list[[1]] = c(seq(x.seg[1,1], x.seg[1,2], 
							x.step), x.seg[1,2])
					}
					y.seg = (y.start + (s$y*height.cn.track))/y.max
					x.seg.circ = c()
					y.seg.circ = c()
					for(k in 1:length(x.seg.list)){
						tmp.x = x.seg.list[[k]]/x.max
						tmp.y = y.seg[k]
						x.seg.circ = c(x.seg.circ , sin(tmp.x*2*pi)*tmp.y)
						y.seg.circ = c(y.seg.circ , cos(tmp.x*2*pi)*tmp.y)
					}
					lines(x.seg.circ, y.seg.circ, col="grey30", lwd=8, 
						lend=1)
				}
				x0 = gen.ref$glob.endpos[nrow(gen.ref)] + gap
				x1 = x0 + indx.gap
				x = c(seq(x0, x1, x.step), x1)
				n.x = length(x)
				x = c(x, x[order(x, decreasing=T)])
				y.start = y.track.start + sum(unlist(n.score.tracks)) + 
					module.gap*n.modules
				y.end = y.start+height.cn.track
				y = c(y.start, y.end)/y.max
				y = rep(y, each=n.x)				
				x.circ = sin((x/x.max)*2*pi)*y
				y.circ = cos((x/x.max)*2*pi)*y
				polygon(x=x.circ, y=y.circ, border="grey50", 
					col= "grey95", lwd=0.5)
				text(0,(y.start+0.7)/y.max, labels="CN", col="grey50", font=2, cex=7)
				# Add chromosome numbers
				nchrom = nrow(gen.ref)
				x = ((gen.ref$glob.startpos[seq(1,nchrom,2)] + 
					gen.ref$glob.endpos[seq(2,nchrom,2)])/2)/x.max
				y = (y.track.start-0.35)/y.max
				rot.vec = -360*x
				rot.vec[rot.vec < -90 & rot.vec > -270] = rot.vec[rot.vec < 
					-90 & rot.vec > -270] + 180
				x.circ = sin((x)*2*pi)*y
				y.circ = cos((x)*2*pi)*y
				n.x = length(x)
				for(j in 1:n.x){
					text(x.circ[j], y.circ[j], labels=chrom[j], col="grey50", 
						font=2, cex=5, srt=rot.vec[j])
				}
				# Make heatmap
				nindx = length(indx.names[[1]])
				indx.col = heat.list
				x.start = 1.1
				x.width = 0.2
				y.start = -0.15
				y.gap = 0.01
				for(j in 1:nindx){
					tmp.col = indx.col[[j]]
					for(k in 1:length(tmp.col)){
						x.line.start = x.start + (k-1)/length(tmp.col)*x.width
						x.line.end = x.line.start + 1.02/
							length(tmp.col)*x.width
						y.line.start = y.start + (j-1)/20
						y.line.end = y.line.start + 0.05 - y.gap
						rect(x.line.start, y.line.start, x.line.end, 
							y.line.end, col=tmp.col[k], border=NA)
					}
					x.lab = x.line.end + 0.01				
					y.lab = (y.line.start + y.line.end)/2-0.007
					text(x.lab, y.lab, labels=indx.names[[1]][j], 
						col="grey40", cex=5, pos=4, font=2)
				}
				y.grid = y.line.end + 0.005
				segments(x.start, y.grid, x.start+ x.width, lwd=8, 
					lend=1)
				x.grid = x.start + x.width*c(0.006, 0.5, 0.994)
				segments(x.grid, y.grid, x.grid, y.grid+0.01, 
					lwd=8, lend=1)
				text(x.grid, rep(y.grid)+0.01, labels=c(0, 0.5, 1), 
					cex=4, pos=3)
				text(0,0, samp.name, cex=10, font=2)
			dev.off()
	}	
}
