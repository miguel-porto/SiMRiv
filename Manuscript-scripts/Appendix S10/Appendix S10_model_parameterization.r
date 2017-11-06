################################################################################
#         EXAMPLES OF INPUT PARAMETER APPROXIMATION FROM REAL DATASETS         #
# This script produces the plots presented in Appendix S2 of the paper         #
#------------------------------------------------------------------------------#
# IMPORTANT NOTE: Run it with the command source('filename'),                  #
# NOT with copy/paste!                                                         #
# In Windows, just drag the file to the R window.                              #
#                                                                              #
################################################################################

library(adehabitatLT)	# for Lévy walk simulation
library(moveHMM)		# for moveHMM simulation
library(TeachingDemos)	# for subplots inside plot

# number of generations to run optimization algorithm
# NOTE: we don't need 1000 generations here, as the algorithm generally converges
# faster. You can reduce this to e.g. 500
generations <- seq(5, 1000, by = 5)
# number of times each candidate solution is simulated to compute average fitness
nrepetitions <- 6
# simulate at 50 times higher frequency than real data
downsample <- 50
# use 7-bin histograms for turning angle and step length ditributions during
# optimization, for solution fitness evaluation
nbins.hist <- c(7, 7)
# use the logarithm of the step lengths for computing the histograms above
log.step.length <- TRUE

################################
# Load/generate some real data #
################################
switch(menu(c(
	"Elk data (Morales et al., 2004)"
	, "moveHMM simulation"
	, "Classical Lévy walk"
	, "Real otter"
	, "SiMRiv simulation"
	, "Repeat last one")), {
	
# Elk dataset
	datasets <- data(package = "moveHMM")$results[, "Item"]
	if(!("elk_data" %in% datasets))
		stop("The elk dataset was not found in package moveHMM. Do you have the latest version of moveHMM installed? If not, type update.packages()")
	data(elk_data)
	filename <- "elk"
# pick coordinates of one elk
	real.data <- elk_data[elk_data[, 1] == "elk-363", 2:3]
# homogeneous landscape
	resistance <- NULL
# so we don't need to specify initial coordinates nor headings
	startCoord <- NULL

# this is just to compute real step lengths
	tmp <- sampleMovement(real.data)
# and then the maximum allowable step length of simulations
	max.step.length <- (max(tmp$stat[, "steplengths"]) / downsample) * 2
# the value above is just a quick calculation to aid the algorithm by providing
# a bounded interval. The algorithm does not need a priori intervals for
# parameters, but obviously it is faster if we provide it some guidance.

# Define a species model to adjust
# For this dataset, we define a very flexible species model with 3 Correlated
# Random Walk states, because we have no idea of how the movement is like.
	species.model <- speciesModel("CRW.CRW.CRW.sl", prob.upp = 0.4
		, steplength = max.step.length)
# This model implies "estimating" 12 parameters:
#    - 3 turning angle concentrations, one for each state, [0, 1]
#    - 3 step lengths, one for each state, [0, given value]
#    - switching probability S1 -> S2, [0, 0.4]
#    - switching probability S1 -> S3, [0, 0.4]
#    - switching probability S2 -> S1, [0, 0.4]
#    - switching probability S2 -> S3, [0, 0.4]
#    - switching probability S3 -> S1, [0, 0.4]
#    - switching probability S3 -> S2, [0, 0.4]

} ,{

# moveHMM simulation (for comments see the Elk)
	filename <- "moveHMM"
# generate a moveHMM trajectory
	stepPar <- c(1, 1, 0.5, 3)		# step distribution parameters
	anglePar <- c(0, 0, 0.001, 3)	# angle distribution parameters 
# specify regression coefficients for the transition probabilities
	beta <- matrix(c(-5, -4), nrow = 1, byrow = TRUE)
# simulate 500 steps of a 2-state movement
	rawData <- as.data.frame(
		simData(nbAnimals = 1, nbStates = 2, stepDist = "gamma", angleDist = "vm",
           stepPar = stepPar, anglePar = anglePar, beta = beta,
           covs = NULL, obsPerAnimal = 500, states = TRUE)
       )

	# extract coordinates
	real.data <- rawData[, 4:5]
	resistance <- NULL
	startCoord <- NULL

	tmp <- sampleMovement(real.data)
	max.step.length <- (max(tmp$stat[, "steplengths"]) / downsample) * 2
	
# a 3-state model better fits in this case (it'll accomodate to 2 states anyway)
	species.model <- speciesModel("CRW.CRW.CRW.sl", prob.upp = 0.4
		, steplength = max.step.length)
	
}, {

# Lévy walk simulation (for comments see the Elk)
	filename <- "levy"
# generate a Lévy walk
	real.data <- simm.levy(1:500)[[1]][, 1:2]
	
	tmp <- sampleMovement(real.data)
	resistance <- NULL
	startCoord <- NULL
	
	max.step.length <- (max(tmp$stat[, "steplengths"]) / downsample) * 2
	
# We don't need more than two states here.
# This 2-state model implies estimating 6 parameters
	species.model <- speciesModel("CRW.CRW.sl", prob.upp = 0.4
		, steplength = max.step.length)
# see attr(species.model, "param.names") to see the parameters to be estimated

}, {

# Load real otter data and resistance raster
# Otters move mainly along rivers, so we have to provide a resistance raster with
# the rivers
# NOTE: because the optimization algorithm is not yet fully optimized for speed,
# running this procedure takes about 22 h! The code is provided just to illustrate
# how to do it.
	if(!file.exists("Appendix S10_otter-realdata.rdata")) {
		stop("Can't find data file.\n************************************************\nPlease copy 'Appendix S10_otter-realdata.rdata' to the folder\n", getwd())
	}
	
	cat("NOTE: this options takes >20 hours to run, see script comments.")
	flush.console()

# let's just downsample 25 times, to make this run faster... but it should work
# equally well with 50	
	downsample <- 25
	
	load("Appendix S10_otter-realdata.rdata")
	filename <- "otter"
	tmp <- sampleMovement(real.data)
	max.step.length <- (max(tmp$stat[, "steplengths"]) / downsample) * 2

# It's a river landscape, so set the starting position to the same as real, to
# be comparable
	startCoord <- matrix(real.data[1, ], nrow = 1)

	species.model <- speciesModel("CRW.CRW.CRW.sl", prob.upp = 0.4
		, steplength = max.step.length, perceptual.range = 50)
# we could specify a simpler model in alternative
#	species.model <- speciesModel("CRW.RW.Rest.sl", prob.upp = 0.4,
# 		steplength = max.step.length, perceptual.range = 50)

}, {

# self simulation - is it able to recover the original simulation parameters
# after 50 times downsampling?
# the answer is "pretty much" - the differences in parameters result in
# negligible differences in the final movement patterns
	filename <- "SiMRiv"
	resistance <- NULL
	startCoord <- NULL

# generate "real" data with SiMRiv
	simsp <- species(state.RW() + state.CRW(0.94), transitionMatrix(0.002, 0.005))
	original <- simulate(simsp, 15000)
	real.data <- sampleMovement(original, downsample)$relocs
	
	tmp <- sampleMovement(real.data)
	max.step.length <- (max(tmp$stat[, "steplengths"]) / downsample) * 2
	
	species.model <- speciesModel("CRW.CRW.sl", prob.upp = 0.4
		, steplength = max.step.length)
	
}, {
	# repeat last one
})

# show the real data for which we are approximating the parameters
plot(real.data, type = "l", asp = 1, main = "Real data")

########################
# Conduct optimization #
########################

time <- system.time({
	solutions <- adjustModel(as.matrix(real.data), species.model
		, resistance = resistance, coords = startCoord
		, resol = downsample, generations = generations
		, nrepetitions = nrepetitions, nbins.hist = nbins.hist
		, step.hist.log = log.step.length, trace = TRUE, TA.variation = FALSE
		, aggregate = TRUE)
})
# ... And that's all it takes to approximate parameters from real data.
# We set trace = TRUE so the values of the objectives for each candidate solution
# are printed every 5 generations.

print(time)		# how long did it take?

# Now visualize results with some fancy plots.

###################################################
# Make optimization plot along generations        #
# good (but not ideal) for assessing convergence  #
###################################################
# Note: for a quick visualization, just type
# generationPlot(solutions, species.model, only.pareto=TRUE)
# the code below is just to beautify the plot for publication!

tiff(paste0("optimization-convergence-", filename, ".tif"), width = 7, height = 7
	, unit = "cm", res = 600, comp = "lzw", pointsize = 9)
par(tcl = -0.25, cex.axis = 0.85)
generationPlot(solutions, species.model, mar = c(5.0, 2.3, 0.9, 2.3)
	, mgp = c(1.2, 0.2, 0), lwd = 0.75, show.legend = FALSE, only.pareto = TRUE)
title(main = "Optimization solutions", font.main = 1)
plot.colors <- rep(c("#000000", "#00bb00", "#ff0000", "#ff9900", "#0077ff"
	, "#9900ff"), 2)
plot.lty <- c(rep(1, 6), rep(2, 6))
if(attr(species.model, "npars") == 6) {
	leg.text <- c("Turn.Ang.conc. S1", "Turn.Ang.conc. S2"
		, expression(paste("Prob. ", S1%->%S2))
		, expression(paste("Prob. ", S2%->%S1))
		, "Step Len. S1", "Step Len. S2")
	legend("bottomleft", xpd = TRUE, inset = c(-0.1, -0.40)
		, legend = leg.text[1:4], lwd = 1, col = plot.colors[1:4], bty = "n"
		, cex = 0.6, lty = plot.lty[1:4])
	legend("bottomright", xpd = TRUE, inset = c(-0.1, -0.40)
		, legend = leg.text[5:6], lwd = 1, col = plot.colors[5:6], bty = "n"
		, cex = 0.6, lty = plot.lty[5:6])
}
if(attr(species.model, "npars") == 12) {
	leg.text <- c("Turn.Ang.conc. S1", "Turn.Ang.conc. S2", "Turn.Ang.conc. S3"
		, expression(paste("Prob. ", S1%->%S2))
		, expression(paste("Prob. ", S1%->%S3))
		, expression(paste("Prob. ", S2%->%S1))
		, expression(paste("Prob. ", S2%->%S3))
		, expression(paste("Prob. ", S3%->%S1))
		, expression(paste("Prob. ", S3%->%S2))
		, "Step Len. S1", "Step Len. S2", "Step Len. S3")
	legend("bottomleft", xpd = TRUE, inset = c(-0.1, -0.40)
		, legend = leg.text[1:4], lwd = 1, col = plot.colors[1:4], bty = "n"
		, cex = 0.6, lty = plot.lty[1:4])
	legend("bottom", xpd = TRUE, inset = c(-0.1, -0.40), legend = leg.text[5:8]
		, lwd = 1, col = plot.colors[5:8], bty = "n", cex = 0.6
		, lty = plot.lty[5:8])
	legend("bottomright", xpd = TRUE, inset = c(-0.1, -0.40)
		, legend = leg.text[9:12], lwd = 1, col = plot.colors[9:12], bty = "n"
		, cex = 0.6, lty = plot.lty[9:12])
}
dev.off()

#######################
# Make movement plots #
#######################
# Just a custom function to plot trajectories and their histograms of turning
# angles and step lengths inset at the corners
movementplot <- function(relocs, downsample = 1, nbins = c(4, 5)
	, range = quantile(relocs.st$stats[, "steplengths"], probs = c(0, 1))
	, resist = NULL) {
	par(mar = c(0.5, 0.5, 2, 0.5), mgp = c(1.2, 0.2, 0), tcl = -0.25)
	relocs.st <- sampleMovement(relocs, downsample)
	if(is.null(resist)) {
		if(downsample == 1)	{
			plot(relocs, type = "l", asp = 1, axes = FALSE)
		} else {
			plot(relocs, type="l", asp=1, col="#aaaaaa", axes = FALSE)
			lines(relocs.st$relocs)
		}
	} else {
		image(resist, axes = FALSE, col = c("white", "#999999"))
		if(downsample == 1)	{
			lines(relocs)
		} else {
			lines(relocs, col = "gray")
			lines(relocs.st$relocs)
		}
	}

	if(nbins[2] > 1) {
		histo <- binCounts(relocs.st$stats[, "steplengths"], range, nbins[2]
			, log = log.step.length)
		subplot(
			barplot(histo, axes = F, main = "Step\nlengths", border = NA
				, names.arg = NA)
			, "bottomright", size = c(0.45, 0.5), type = "plt"
			, pars = list(mar = c(0, 0, 4, 0), cex = 0.7, font.main = 1)
		)
	}

	if(nbins[1] > 1) {
		histota <- binCounts(relocs.st$stats[, "turningangles"], c(-pi, pi)
			, nbins[1], FALSE)
		subplot(
			barplot(histota, axes = FALSE, main = "Turning\nangles", border = NA
				, names.arg = NA)
			, "bottomleft", size = c(0.45, 0.5), type = "plt"
				, pars = list(mar = c(0, 0, 4, 0), cex = 0.7, font.main = 1)
		)
	}
	return(range)
}


# Make 3 example simulations with the optimized parameters
sim <- lapply(1:3, function(i) simulate(species.model(
	solutions[[length(solutions)]]$par[i, ]), nrow(real.data) * downsample
	, resist = resistance, coords = startCoord))

# do plots with real and simulated trajectories, with histograms,
# for visual comparison
tiff(paste0("real-vs-simulated-", filename, ".tif"), width = 8 * 2, height = 8 * 2
	, unit = "cm", res = 600, comp = "lzw", pointsize = 9)
par(mfrow = c(2,2), font.main = 1, cex.main = 1.3)

if(exists("original")) {
	rh <- movementplot(original, downsample, nbins = nbins.hist
		, resist = resistance)
	rm(original)
} else {
	rh <- movementplot(real.data, 1, nbins = nbins.hist, resist = resistance)
}
title(main = "Real")
box(bty = "o", lwd = 0.5)
sapply(sim, function(m) {
	movementplot(m, downsample, range = rh, nbins = nbins.hist, resistance)
	title(main = "Simulated")
	box(bty = "o", lwd = 0.5)
})
dev.off()

