plot.relocs = function(relocs, resist, n.rep = NA , xlim = NA, ylim = NA
	, kde, col = ifelse(is.na(n.rep), "#0000ff", "#0000ff30"), lwd = 0.2, scale = F, contour.lines = c(5, 25, 50, 75, 95)) {
	if(all(is.na(xlim))) xlim = range(relocs[, 1])
	if(all(is.na(ylim))) ylim = range(relocs[, 2])

	if(!missing(resist) && !is.null(resist) && !is.na(resist)) {
		image(resist, maxpixels=Inf, asp=1, zlim=c(0,1), col=c(gray(seq(1, 0.6, len=20)), "#990000")
			, ylim = ylim, xlim = xlim, axes = F, xlab = NA, ylab = NA)
	} else {
		plot.new()
		plot.window(xlim = xlim, ylim = ylim, asp = 1)
	}
	if(is.na(n.rep)) {
		lines(relocs[, 1:2], lwd = lwd, col = col)
#		segments(relocs[from, 1], relocs[from, 2], relocs[to, 1], relocs[to, 2]
#			, lwd = 0.2, col = col) #col=c("#ff0000","#0000ff")[relocs[,3]+1],lwd=1)
	} else if(n.rep == 1) {
		from = 1:(nsteps - 1)
		to = from + 1
		segments(relocs[from, 1], relocs[from, 2], relocs[to, 1], relocs[to, 2]
			, lwd = lwd, col = col)
	} else {
		nsteps = nrow(relocs) / n.rep
		for(burst in 0:(n.rep-1)) {
			subset = (burst * nsteps + 1):((burst + 1) * nsteps)
			lines(relocs[subset, 1:2], lwd = lwd, col = col)
#			segments(relocs[from, 1], relocs[from, 2], relocs[to, 1], relocs[to, 2]
#				, lwd = 0.2, col = col)
		}
	}
#plot(kde1, disp="slice", cont=c(80, 25, 5), col=c("#00bb00", "#00dd00", "#00ff00"), add=T, lwd=7, drawlabels=F)
	if(!missing(kde) && !is.na(kde)) {
		n = length(contour.lines)
		plot(kde, disp="slice", cont = contour.lines, col = "#00000044", add = T, lwd = 2, drawlabels = F)
		#plot(kde, disp="slice", cont = contour.lines, col = rgb(((0:(n - 1)) + 0) / (n - 1), ((n:1) + n) / (n * 2), 0), add = T, lwd = 0.75, drawlabels = F)
		plot(kde, disp="slice", cont = contour.lines, col = c("#ff6600", "#ffff00"), add = T, lwd = 0.75, drawlabels = F)
	}
	
	if(scale) {
		usr = par("usr")
		xw = usr[2] - usr[1]
		yw = usr[4] - usr[3]
		x0 = usr[1] + xw * 0.05
		y0 = usr[3] + yw * 0.05
		x1 = x0 + 1000
		segments(x0, y0, x1, y0, lwd = 5, lend = "butt")
		text((x0 + x1) / 2, y0 + yw * 0.02, "1 km", adj = c(0.5, 0))
	}
	
	box(bty = "o")
}


