#################################################################################
#  This script produces the plots of figures 2-4 presented in the main text     #
# of the paper                                                                  #
#-------------------------------------------------------------------------------#
# IMPORTANT NOTE: Run it with the command source('filename'), v                 #
# NOT with copy/paste! In Windows, just drag the file to the R window.			#
#                                                                               #
# This script must have access to the data file 'resistance-rasters.rdata'.     #
# Either place this rdata file in the default working directory, or use setwd() #
# to set a working directory.                                                   #
#                                                                               #
# Figs. 3 & 4 take about 10 minutes to run in 8 cores                           #
#################################################################################

library(SiMRiv)
library(ks)			# 2D kernel density
library(parallel)	# for parallel lapply

# nr of steps to simulate (10000 in the paper)
nsteps <- 10000
# nr of repeated simulations to conduct in each analysis (1000 in the paper)
n.rep <- 1000
# the step length to be used in simulations with the river raster
step.length <- 10
# number of cores to use in parallel simulations; adjust as necessary
n.cores <- 8

if(!file.exists("Appendix S9_resistance-rasters.rdata")) {
	stop("Can't find data file.\n************************************************\nPlease copy 'Appendix S9_resistance-rasters.rdata' to the folder\n", getwd())
}

##  Load rasters corresponding to how different animals might see the landscape:
# - terrestrial animal: Moves in any land use except urban and avoids water
#   (but can swim). e.g. wolf
# - semi-aquatic animal: Moves primarily in water, but may move overland at times,
#   more often in forested areas. e.g. otter
# - aquatic animal: Moves exclusively in water. e.g. fish

load("Appendix S9_resistance-rasters.rdata")
# this is a list with 3 elements: terr, amph, fish, corresponding to the
# descriptions above

# A custom wrapper for repeated parallelized simulations
my.simulate <- function(resist, nrep, species, init = NULL) {
	stack <- matrix(ncol = 3, nrow = 0)
	print(species); flush.console()
	cl <- makeCluster(n.cores)
	clusterExport(cl, c("xyFromCell", "values", "resist", "species", "nsteps"
		, "simulate"), envir = environment())
	stack <- parLapply(cl, 1:nrep, function(i) {
		if(is.null(init)) {
			if(is.null(resist)) # start at (0,0)
				init <- c(0, 0)
			else # start in a random point with zero resistance
				init <- xyFromCell(resist, sample(which(values(resist) == 0), 1))
		}
		relocs <- simulate(species, nsteps, resist = resist, coords = init)
		return(relocs)
	})
	stopCluster(cl)
# join results in one matrix
	stack <- do.call(rbind, stack)
	cat("\n")
	return(stack)
}

# A custom function to plot movements and resistance raster styled like the
# figures in the main text
plot.relocs <- function(resist, relocs, n.rep = NA , xlim = NA, ylim = NA
	, kde, col = "#0000ff30", lwd = 0.2, scale = FALSE
	, contour.lines = c(5, 25, 50, 75, 95)) {
	
	if(all(is.na(xlim))) xlim <- range(relocs[, 1])
	if(all(is.na(ylim))) ylim <- range(relocs[, 2])

	if(!is.null(resist)) {
		image(resist, maxpixels = Inf, asp = 1, zlim = c(0, 1)
			, col = c(gray(seq(1, 0.6, len = 20)), "#990000"), ylim = ylim
			, xlim = xlim, axes = FALSE, xlab = NA, ylab = NA)
	} else {
		plot.new()
		plot.window(xlim = xlim, ylim = ylim, asp = 1)
	}
	if(is.na(n.rep)) {
		lines(relocs[, 1:2], lwd = lwd, col = col)
	} else if(n.rep == 1) {
		from <- 1:(nsteps - 1)
		to <- from + 1
		segments(relocs[from, 1], relocs[from, 2], relocs[to, 1], relocs[to, 2]
			, lwd = lwd, col = col)
	} else {
		nsteps <- nrow(relocs) / n.rep
		for(burst in 0:(n.rep - 1)) {
			subset <- (burst * nsteps + 1):((burst + 1) * nsteps)
			lines(relocs[subset, 1:2], lwd = lwd, col = col)
		}
	}

	if(!missing(kde) && !is.na(kde)) {
		n <- length(contour.lines)
		plot(kde, disp="slice", cont = contour.lines, col = "#00000044"
			, add = TRUE, lwd = 2, drawlabels = FALSE)
		plot(kde, disp="slice", cont = contour.lines, col = c("#ff6600", "#ffff00")
			, add = TRUE, lwd = 0.75, drawlabels = FALSE)
	}
	
	if(scale) {
		usr <- par("usr")
		xw <- usr[2] - usr[1]
		yw <- usr[4] - usr[3]
		x0 <- usr[1] + xw * 0.05
		y0 <- usr[3] + yw * 0.05
		x1 <- x0 + 1000
		segments(x0, y0, x1, y0, lwd = 5, lend = "butt")
		text((x0 + x1) / 2, y0 + yw * 0.02, "1 km", adj = c(0.5, 0))
	}
	
	box(bty = "o")
}

## START MENU
choice <- menu(c(
	"Fig.3 - Influence of landscape"
	, "Fig.4 - Influence of perceptual range"
	, "Fig.2 - Basic plots"
	, "Just plot last run"
))

switch(choice , {
	output.filename <- "Figure-3-influence-landscape"
	simulated.rasters <- rasters
	# create species
	levy.walker <- species((state.RW() * 100) + (state.CRW(0.95) * 500)
		, transitionMatrix(0.01, 0.002)) + step.length

	cat("Simulating...\n"); flush.console()
	
	time <- system.time({
		relocs <- mapply(my.simulate, simulated.rasters, n.rep, list(levy.walker)
			, SIMPLIFY = FALSE)
	})
	
	cat("Took", time["elapsed"], "secs.\n")
		
	config <- list(
		list(title = "Terrestrial")
		, list(title = "Semi-aquatic")
		, list(title = "Aquatic")
	)
	
	# a complex layout for the multi-panel figure
	layout <- matrix(c(8,1,1,9,2,2,10,3,3,8,4,5,9,6,6,10,7,7), ncol = 2
		, byrow = FALSE)
	heights <- c(1,4,4,1,4,4,1,4,4)
	
	# these are the bounding boxes of the close-ups in the figure
	xlim <- list(c(688558, 695000), c(688558, 695000), c(688558, 695000)
			, c(688577.7, 690937.3), c(690432.9, 694948.1), c(692000, 693600)
			, c(689100, 690600))
	ylim <- list(c(4622621, 4628370), c(4622621, 4628370), c(4622621, 4628370)
				, c(4626924.0, 4627964.0), c(4622351.0, 4624341.0)
				, c(4624289, 4625711), c(4626834, 4628166))
	raster.indices <- c(1,2,3,1,1,2,3)
	cols <- c(rep("#0000ff20", 3), rep("#0000ff40", 4))
	corners <- list("A)", "B)", "C)", expression(paste(Ai, ")"))
		, expression(paste(Aii, ")")), expression(paste(Bi, ")"))
		, expression(paste(Ci, ")")))
	
	rects <- list(
		list(c(xlim[[4]], ylim[[4]]), c(xlim[[5]], ylim[[5]]))
		, list(c(xlim[[6]], ylim[[6]]))
		, list(c(xlim[[7]], ylim[[7]]))
		, NULL, NULL, NULL, NULL
	)

}, {
	output.filename <- "Figure-4-influence-percep-range"
	species <- lapply(c(5000, 2000, 500), function(pr) {
		species(list(
			state(0, perceptualRange("cir", pr), step.length)
			,state(0.95, perceptualRange("cir", pr), step.length)
		), transitionMatrix(0.01, 0.002), paste("P.Range", pr))
	})
	
	simulated.rasters <- rasters[rep("amph", 3)]
	cat("Simulating...\n"); flush.console()

	time <- system.time({
		relocs <- mapply(my.simulate, simulated.rasters, n.rep, species
			, SIMPLIFY = FALSE)
	})
	
	cat("Took", time["elapsed"], "secs.\n")
	
	config <- list(
		list(title = "5000 m")
		, list(title = "2000 m")
		, list(title = "500 m")
	)
	layout <- matrix(c(7,1,8,2,9,3,7,4,8,5,9,6), ncol = 2, byrow = FALSE)
	heights <- c(1,8,1,8,1,8)
	
	xlim <- list(c(688558, 695000), c(688558, 695000), c(688558, 695000)
			, c(688558, 695000), c(688558, 695000), c(688558, 695000))
	ylim <- list(c(4622621, 4628370), c(4622621, 4628370), c(4622621, 4628370)
			, c(4622621, 4628370), c(4622621, 4628370), c(4622621, 4628370))

	xlim <- list(c(687241, 699253) - 500, c(687241, 699253) - 500, c(687241, 699253) - 500
			, c(692428, 695547) - 500, c(692428, 695547) - 500, c(692428, 695547) - 500)
	ylim <- list(c(4630366, 4640317), c(4630366, 4640317), c(4630366, 4640317)
			, c(4634402, 4636986), c(4634402, 4636986), c(4634402, 4636986))

	raster.indices <- c(1:3, 1:3)
	cols <- c(rep("#0000ff20", 3), rep("#0000ff40", 3))
	corners <- list("A)", "B)", "C)", expression(paste(Ai, ")"))
		, expression(paste(Bi, ")")), expression(paste(Ci, ")")))
	rects <- list(
		list(c(xlim[[4]], ylim[[4]]))
		, list(c(xlim[[5]], ylim[[5]]))
		, list(c(xlim[[6]], ylim[[6]]))
		, NULL, NULL, NULL
	)

}, {	# comparison river / homogeneous
	n.rep <- 1
	nsteps <- 3000
	simulated.rasters <- list(NULL, rasters[["amph"]], NULL, rasters[["amph"]])
	raster.indices <- c(NA,2,NA,2)
	pR <- 200	# perceptual range

	species <- list(	
		species(state.RW()) * pR + step.length
		, species(state.RW()) * pR + step.length
		, species(state.RW() + state.CRW(0.95)
			, transitionMatrix(0.01, 0.01)) * pR + step.length
		, species(state.RW() + state.CRW(0.95)
			, transitionMatrix(0.01, 0.01)) * pR + step.length
	)
	config <- list(
		list(title = "Random walker\nhomogeneous landscape")
		, list(title = "Random walker\nriver")
		, list(title = "2-state Lévy-like walker\nhomogeneous landscape")
		, list(title = "2-state Lévy-like walker\nriver")
	)
	
	relocs <- mapply(my.simulate, simulated.rasters, n.rep, species
		, list(c(691391, 4629949)), SIMPLIFY = FALSE)
	
	tiff("Figure-2-basic-illustration.tif", width = 10 * 2, height = 10 * 2
		, unit = "cm", comp = "lzw", res = 1000)
	par(mar = c(0.1, 0.1, 2.5, 0.1))
	layout(matrix(c(seq_along(relocs)), nrow = 2, byrow = TRUE))
	
	for(i in seq_along(relocs)) {
		if(n.rep == 1)
			col <- relocs[[i]][, 3] + 1
		else
			col <- "#0000ff20"
		plot.relocs(simulated.rasters[[i]], relocs[[i]], n.rep
			, xlim = NA, ylim = NA, col = col, lwd = 1, scale = FALSE)
		title(main = config[[i]]$title, font.main = 1)
		mtext(paste0(LETTERS[[i]], ")"), line = -1.8, adj = 0.025, padj = 0, cex = 1.3)
	}
	dev.off()
}, {
	# just plot last run
})


## Make Fig.3 & 4, and the full res plots in appendix
if(choice %in% c(1, 2)) {
	cat("Computing kernel densities...\n"); flush.console()
	if(n.rep > 1) {
		kde <- lapply(relocs, function(rel) {
			kde(rel[, 1:2], H = matrix(c(15000, 0, 0, 15000), nc = 2)
				, binned = TRUE, bgridsize = c(1000,1000))
		})
	} else {
		kde <- lapply(relocs, function(x) NA)
	}

	cat("Making plots...\n"); flush.console()

	## Combined excerpt plot
	tiff(paste0(output.filename, ".tif"), width = 10 * 2, height = 10 * 3
		, unit = "cm", comp = "lzw", res = 400)
	par(mar = c(0.1, 0.1, 0.1, 0.1))
	layout(layout, heights = heights)
	for(i in seq_along(raster.indices)) {
		plot.relocs(simulated.rasters[[raster.indices[i]]]
			, relocs[[raster.indices[i]]], n.rep, xlim = xlim[[i]]
			, ylim = ylim[[i]], col = cols[i], lwd = 0.75, scale = TRUE)
		mtext(corners[[i]], line = -1.8, adj = 0.025, padj = 0, cex = 1.3)
		if(!is.null(rects[[i]])) {
			for(j in seq_along(rects[[i]])) {
				r <- rects[[i]][[j]]
				rect(r[1], r[3], r[2], r[4], border = "black", lty = 2)
			}
		}
	}
	# now the labels
	par(mar = rep(0, 4))
	for(i in seq_along(config)) {
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(0.5, 0.5, config[[i]]$title, adj = c(0.5, 0.5), cex = 2)
	}
	dev.off()

	## Full extent high-res plots
	for(i in seq_along(relocs)) {
		tiff(paste0(output.filename, "-full-", i, ".tif"), width = 21
			, height = 21 * diff(range(relocs[[i]][, 2])) / diff(range(relocs[[i]][, 1]))
			, unit = "cm", comp = "lzw", res = 800)
		par(mar = c(0.1, 0.1, 1.5, 0.1))
		if(n.rep == 1)
			col <- relocs[[i]][, 3] + 1
		else
			col <- "#0000ff20"
		plot.relocs(simulated.rasters[[i]], relocs[[i]], n.rep, kde = kde[[i]]
			, col = col, lwd = 0.2, scale = TRUE, contour.lines = c(25, 95))
		title(main = config[[i]]$title)
		dev.off()
	}
}

