## Define a species model to fit to real data
# Species models differ in the number of states and in which parameters are
# allowed to vary
speciesModel <- function(type, perceptual.range = 0, steplength = 1
	, prob.upperbound = 0.5, max.concentration = 0.99) {
	
	return(switch(pmatch(type, c("CRW", "CRW.sl", "RW.CRW", "CRW.CRW", "CRW.pw"
		, "RW.CRW.sl", "CRW.CRW.sl", "CRW.CRW.CRW.sl", "CRW.CRW.CRW.sl.symm", "CRW.RW.Rest.sl")), {
		f <- function(parameters) {
			return(species(
				state.CRW(parameters[1])
			) * perceptual.range + steplength)
		}
		attr(f, "npars") <- 1
		attr(f, "lower.bounds") <- 0
		attr(f, "upper.bounds") <- max.concentration
		attr(f, "param.names") <- "Turning angle concentration"
		attr(f, "param.types") <- c("TA1")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				state.CRW(parameters[1])
			) * perceptual.range + parameters[2])
		}
		attr(f, "npars") <- 2
		attr(f, "lower.bounds") <- rep(0, 2)
		attr(f, "upper.bounds") <- c(max.concentration, steplength)
		attr(f, "param.names") <- c("Turning angle concentration", "Max step length")
		attr(f, "param.types") <- c("TA1", "SL1")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				state.CRW(parameters[1]) + state.RW()
				, transitionMatrix(parameters[2], parameters[3])
			) * perceptual.range + steplength)
		}
		attr(f, "npars") <- 3
		attr(f, "lower.bounds") <- c(0, rep(0, 2))
		attr(f, "upper.bounds") <- c(max.concentration, rep(prob.upperbound, 2))
		attr(f, "param.names") <- c("Turning angle concentration", "Prob. CRW -> RW"
			, "Prob. RW -> CRW")
		attr(f, "param.types") <- c("TA1", "12", "21")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				state.CRW(parameters[1]) + state.CRW(parameters[2])
				, transitionMatrix(parameters[3], parameters[4])
			) * perceptual.range + steplength)
		}
		attr(f, "npars") <- 4
		attr(f, "lower.bounds") <- c(0, 0, rep(0, 2))
		attr(f, "upper.bounds") <- c(max.concentration, max.concentration, rep(prob.upperbound, 2))
		attr(f, "param.names") <- c("Turning angle concentration S1"
			, "Turning angle concentration S2", "Prob. S1 -> S2"
			, "Prob. S2 -> S1")
		attr(f, "param.types") <- c("TA1", "TA2", "12", "21")
		return(f)
	}, {
		if(perceptual.range <= 1)
			stop("You must provide the maximum allowable perceptual range size, e.g. perceptual.range = 500")
		f <- function(parameters) {
			return(species(
				state.CRW(parameters[1])
			) * parameters[2] + steplength)
		}
		attr(f, "npars") <- 2
		attr(f, "lower.bounds") <- c(0, 1)
		attr(f, "upper.bounds") <- c(max.concentration, perceptual.range)
		attr(f, "param.names") <- c("Turning angle concentration"
			, "Perceptual range radius")
		attr(f, "param.types") <- c("TA1", "PR1")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				(state.CRW(parameters[1]) + parameters[5]) + (state.RW() + parameters[4])
				, transitionMatrix(parameters[2], parameters[3])
			) * perceptual.range)
		}
		attr(f, "npars") <- 5
		attr(f, "lower.bounds") <- c(0, rep(0, 4))
		attr(f, "upper.bounds") <- c(max.concentration, rep(prob.upperbound, 2), rep(steplength, 2))
		attr(f, "param.names") <- c("Turning angle concentration"
			, "Prob. CRW -> RW", "Prob. RW -> CRW", "Max step length RW"
			, "Max step length CRW")
		attr(f, "param.types") <- c("TA1", "12", "21", "SL2", "SL1")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				(state.CRW(parameters[1]) + parameters[5]) 
					+ (state.CRW(parameters[2]) + parameters[6])
				, transitionMatrix(parameters[3], parameters[4])
			) * perceptual.range)
		}
		attr(f, "npars") <- 6
		attr(f, "lower.bounds") <- rep(0, 6)
		attr(f, "upper.bounds") <- c(max.concentration, max.concentration
			, rep(prob.upperbound, 2), rep(steplength, 2))
		attr(f, "param.names") <- c("Turning angle concentration S1"
			, "Turning angle concentration S2", "Prob. S1 -> S2", "Prob. S2 -> S1"
			, "Max step length S1", "Max step length S2")
		attr(f, "param.types") <- c("TA1", "TA2", "12", "21", "SL1", "SL2")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				(state.CRW(parameters[1]) + parameters[10])
					+ (state.CRW(parameters[2]) + parameters[11])
					+ (state.CRW(parameters[3]) + parameters[12])
				, transitionMatrix(parameters[4], parameters[5], parameters[6]
					, parameters[7], parameters[8], parameters[9])
			) * perceptual.range)
		}
		attr(f, "npars") <- 12
		attr(f, "lower.bounds") <- rep(0, 12)
		attr(f, "upper.bounds") <- c(rep(max.concentration, 3)
			, rep(prob.upperbound, 6), rep(steplength, 3))
		attr(f, "param.names") <- c("Turning angle concentration S1"
			, "Turning angle concentration S2", "Turning angle concentration S3"
			, "Prob. S1 -> S2", "Prob. S1 -> S3", "Prob. S2 -> S1"
			, "Prob. S2 -> S3", "Prob. S3 -> S1", "Prob. S3 -> S2"
			, "Max step length S1", "Max step length S2", "Max step length S3")
		attr(f, "param.types") <- c("TA1", "TA2", "TA3", "12", "13", "21", "23"
			, "31", "32", "SL1", "SL2", "SL3")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				(state.CRW(parameters[1]) + parameters[7])
					+ (state.CRW(parameters[2]) + parameters[8])
					+ (state.CRW(parameters[3]) + parameters[9])
				, transitionMatrix(parameters[4], parameters[5], parameters[4]
					, parameters[6], parameters[5], parameters[6])
			) * perceptual.range)
		}
		attr(f, "npars") <- 9
		attr(f, "lower.bounds") <- rep(0, 9)
		attr(f, "upper.bounds") <- c(rep(max.concentration, 3)
			, rep(prob.upperbound, 3), rep(steplength, 3))
		attr(f, "param.names") <- c("Turning angle concentration S1"
			, "Turning angle concentration S2", "Turning angle concentration S3"
			, "Prob. S1 <-> S2", "Prob. S1 <-> S3", "Prob. S2 <-> S3"
			, "Max step length S1", "Max step length S2", "Max step length S3")
		attr(f, "param.types") <- c("TA1", "TA2", "TA3", "12", "13", "23"
			, "SL1", "SL2", "SL3")
		return(f)
	}, {
		f <- function(parameters) {
			return(species(
				(state.CRW(parameters[1]) + parameters[6])
					+ (state.RW() + parameters[7]) + state.Resting()
				, transitionMatrix(parameters[2], 0, parameters[3], parameters[4]
					, 0, parameters[5])
			) * perceptual.range)
		}
		attr(f, "npars") <- 7
		attr(f, "lower.bounds") <- rep(0, 7)
		attr(f, "upper.bounds") <- c(max.concentration, rep(prob.upperbound, 4)
			, rep(steplength, 2))
		attr(f, "param.names") <- c("Turning angle concentration CRW"
			, "Prob. CRW -> RW", "Prob. RW -> CRW", "Prob. RW -> Rest"
			, "Prob. Rest -> RW", "Max step length CRW", "Max step length RW")
		attr(f, "param.types") <- c("TA1", "12", "21", "23", "32", "SL1", "SL2")
		return(f)
	}))
}

## Find numerical approximations of input paramters for a read given trajectory
# and a species model, at the given temporal frequency.
# See manual for explanation of the function parameters
adjustModel <- function(
	realData
	, species.model
	, resolution = 10
# simulation parameters
	, resistance = NULL
	, coords = NULL
	, angles = NULL
# fitness function parameters
	, nrepetitions = 6
	, nbins.hist = if(aggregate.obj.hist) c(7, 7, 0) else c(3, 3, 0)
	, step.hist.log = TRUE
	, nlags = 100
	, window.size = dim(reference$stats)[1] %/% nlags
	, aggregate.obj.hist = TRUE
	, step.hist.range = c(0, 1)
# GA options
	, popsize = 100, generations = seq(5, 1000, by=5), mprob = 0.2
# TODO: if using a raster, parallel performance doesn't increase because
# it is loaded in all R processes. This should be fixed.
	, parallel = is.null(resistance)
	, trace = TRUE
) {
	realData <- as.matrix(realData)
	nsteps <- dim(realData)[1]
# compute turning angles and step lengths of the real movement
	reference <- sampleMovement(realData, resist = resistance)
	
	if(nbins.hist[3] > 0) {
# in this case, use the variation in turning angles along a moving window as
# objectives.
# So, compute the SD of turning angles in a fixed-size juxtaposed moving window
		a.var.ref <- angle.variation(reference, window.size = window.size)

# make the (fixed-range) histogram of all the moving window SDs
# the range of the histogram is 10% larger than observed range to account for
# stochasticity in simulations
		increase <- (diff(range(a.var.ref)) * 0.1) / 2
		range.varta <- range(a.var.ref) + c(-increase, increase)
		hist.tv.ref <- binCounts(a.var.ref, range.varta, nbins.hist[3])
	} else {
		range.varta <- NA
		hist.tv.ref <- NA
	}
# otherwise just use the plain old turning angle histogram - trials show that
# this works fine for most cases, so better not to complicate much
	hist.ta.ref <- binCounts(reference$stats[, "turningangles"], c(-pi, pi)
		, nbins.hist[1], FALSE)
	
# make the step length histogram	
	range.step <- quantile(reference$stats[, "steplengths"]
		, probs = step.hist.range)
	hist.sl.ref <- binCounts(reference$stats[, "steplengths"], range.step
		, nbins.hist[2], step.hist.log)

	cl <- NULL
	cat("Real data:", nsteps,", Simulating", nsteps * resolution,"steps @ res"
		, resolution,"\n")
#	cat("Time window size for computing turning angle variation:", window.size, "\n")
	
	if(inherits(parallel, "cluster")) {
		cl <- parallel
	} else if(inherits(parallel, "logical")) {
		if(parallel) {
			cl <- makeCluster(detectCores() - 1)
		}
	} else if(inherits(parallel, "numeric")) {
		cl <- makeCluster(parallel)
	}

	if(!is.null(cl)) {	# GO PARALLEL
		clusterCall(cl, function() library(SiMRiv))
		clusterExport(cl, c("nrepetitions", "nbins.hist", "nsteps"
			, "range.varta", "window.size", "range.step", "resolution"
			, "species.model", "resistance", "coords", "angles"
			), envir = environment())
	
		objective.function.par <- function(inp.mat, reference) {
			crit <- parApply(cl, inp.mat, 1, function(inp.par, ref) {
				sp.sim <- species.model(inp.par)
				
				if(nrepetitions == 1) {
					rel <- simulate(sp.sim, nsteps * resolution
						, resist = resistance, coords = coords, angles = angles)
					s <- sampleMovement(rel, resolution, resist = resistance)
					
					if(nbins.hist[3] > 0) {
						a.var.sim <- angle.variation(s, window.size = window.size)
						hist.tv <- binCounts(a.var.sim, range.varta, nbins.hist[3])
					} else
						hist.tv <- NA
					hist.ta <- binCounts(s$stats[, "turningangles"]
						, c(-pi, pi), nbins.hist[1], FALSE)
					hist.sl <- binCounts(s$stats[, "steplengths"], range.step
						, nbins.hist[2], step.hist.log)
				} else {
					hist.ta.mat <- matrix(ncol = nbins.hist[1], nrow = nrepetitions)
					hist.sl.mat <- matrix(ncol = nbins.hist[2], nrow = nrepetitions)
					hist.tv.mat <- matrix(ncol = nbins.hist[3], nrow = nrepetitions)
					for(i in seq_len(nrepetitions)) {
						rel <- simulate(sp.sim, nsteps * resolution
							, resist = resistance, coords = coords, angles = angles)
						s <- sampleMovement(rel, resolution, resist = resistance)
						if(nbins.hist[3] > 0) {
							a.var.sim <- angle.variation(s
								, window.size = window.size)
							hist.tv.mat[i, ] <- binCounts(a.var.sim, range.varta
								, nbins.hist[3])
						}
						hist.ta.mat[i, ] <- binCounts(s$stats[, "turningangles"]
							, c(-pi, pi), nbins.hist[1], FALSE)
						hist.sl.mat[i, ] <- binCounts(s$stats[, "steplengths"]
							, range.step, nbins.hist[2], step.hist.log)
					}
					hist.ta <- apply(hist.ta.mat, 2, mean)
					hist.sl <- apply(hist.sl.mat, 2, mean)
					hist.tv <- apply(hist.tv.mat, 2, mean)
				}
				
				if(nbins.hist[1] > 0)
					crit.ta <- abs(log(hist.ta + 1) - log(ref[["hist.ta"]] + 1))
				else
					crit.ta <- NULL
				
				if(nbins.hist[2] > 0) {
# We put a log here. This does make a difference when objectives are aggregated
# because a difference of 1 between 2 and 3 should be given more weight than a
# difference of 1 between 20 and 21. If objectives are not aggregated, to log or
# not to log is irrelevant.
# Bins with small values are usually determinant in the final movement pattern.
					crit.sl <- abs(log(hist.sl + 1) - log(ref[["hist.sl"]] + 1))
				} else
					crit.sl <- NULL

				if(nbins.hist[3] > 0)
					crit.tv <- abs(log(hist.tv + 1) - log(ref[["hist.tv"]] + 1))
				else
					crit.tv <- NULL

				if(aggregate.obj.hist) {
# in this case we sum the absolute differences in each histogram bar to use as
# objectives
					crit <- c("TA.diff" = if(is.null(crit.ta)) NULL else sum(crit.ta)
						, "SL.diff" = if(is.null(crit.sl)) NULL else sum(crit.sl)
						, "TV.diff" = if(is.null(crit.tv)) NULL else sum(crit.tv)
					)
				} else
# otherwise we use the bar-wise absolute differences as objectives
					crit <- c(crit.ta, crit.sl, crit.tv)

				return(crit)
			}, reference)
			
			if(is.null(dim(crit))) {
				dim(crit) <- c(1, length(crit))
			}
			
			if(trace) print(round(crit, 1))
			return(crit)
		}
		
# Run NSGA-II optimization
		sol <- nsga2(objective.function.par, attr(species.model, "npars")
			, ifelse(aggregate.obj.hist, sum(nbins.hist > 0), sum(nbins.hist))
# the histograms of real data to compare simulations with
			, list(hist.ta = hist.ta.ref, hist.sl = hist.sl.ref, hist.tv = hist.tv.ref)
			, generations = generations, popsize = popsize
			, lower.bounds = attr(species.model, "lower.bounds")
			, upper.bounds = attr(species.model, "upper.bounds")
			, vectorized = TRUE, mprob = mprob
		)
	} else {	# GO SINGLE CORE
# see comments on the objective function above
		objective.function = function(inp.par, ref) {
			if(trace) {
				cat("INP: ")
				for(o in inp.par)
					cat(sprintf("%5.2f ", o))
				cat("\n")
				flush.console()
			}

			sp.sim = species.model(inp.par)

			if(nrepetitions == 1) {
				rel = simulate(sp.sim, nsteps * resolution, resist = resistance
					, coords = coords, angles = angles)
				s = sampleMovement(rel, resolution, resist = resistance)
				if(nbins.hist[3] > 0) {
					a.var.sim <- angle.variation(s, window.size = window.size)
					hist.tv <- binCounts(a.var.sim, range.varta, nbins.hist[3])
				} else
					hist.tv <- NA
				hist.ta <- binCounts(s$stats[, "turningangles"], c(-pi, pi)
					, nbins.hist[1], FALSE)

				hist.sl <- binCounts(s$stats[, "steplengths"], range.step
					, nbins.hist[2], step.hist.log)
			} else {
				hist.ta.mat <- matrix(ncol = nbins.hist[1], nrow = nrepetitions)
				hist.sl.mat <- matrix(ncol = nbins.hist[2], nrow = nrepetitions)
				hist.tv.mat <- matrix(ncol = nbins.hist[3], nrow = nrepetitions)
		
				for(i in seq_len(nrepetitions)) {
					rel <- simulate(sp.sim, nsteps * resolution
						, resist = resistance, coords = coords, angles = angles)
					s <- sampleMovement(rel, resolution, resist = resistance)
					if(nbins.hist[3] > 0) {
						a.var.sim <- angle.variation(s, window.size = window.size)
						hist.tv.mat[i, ] <- binCounts(a.var.sim, range.varta
							, nbins.hist[3])
					}
					hist.ta.mat[i, ] <- binCounts(s$stats[, "turningangles"]
						, c(-pi, pi), nbins.hist[1], FALSE)

					hist.sl.mat[i, ] <- binCounts(s$stats[, "steplengths"]
						, range.step, nbins.hist[2], step.hist.log)
				}
				hist.ta <- apply(hist.ta.mat, 2, mean)
				hist.sl <- apply(hist.sl.mat, 2, mean)
				hist.tv <- apply(hist.tv.mat, 2, mean)
			}

			if(nbins.hist[1] > 0)
				crit.ta <- abs(log(hist.ta + 1) - log(ref[["hist.ta"]] + 1))
			else
				crit.ta <- NULL
			
			if(nbins.hist[2] > 0) {		# see comment on the log above
				crit.sl <- abs(log(hist.sl + 1) - log(ref[["hist.sl"]] + 1))
			} else
				crit.sl <- NULL

			if(nbins.hist[3] > 0)
				crit.tv <- abs(log(hist.tv + 1) - log(ref[["hist.tv"]] + 1))
			else
				crit.tv <- NULL

			if(aggregate.obj.hist) {
# in this case we sum the absolute differences in each histogram bar to use
# as objectives
					crit <- c("TA.diff" = if(is.null(crit.ta)) NULL else sum(crit.ta)
						, "SL.diff" = if(is.null(crit.sl)) NULL else sum(crit.sl)
						, "TV.diff" = if(is.null(crit.tv)) NULL else sum(crit.tv)
					)
			} else
# otherwise we use the bar-wise absolute differences as objectives
				crit <- c(crit.ta, crit.sl, crit.tv)

			if(trace) {
				cat("OBJ: ")
				for(o in crit)
					cat(sprintf("%5.2f ", o))
				cat("\n")
				flush.console()
			}
			return(crit)
		}

# Run NSGA-II optimization
		sol <- nsga2(objective.function, attr(species.model, "npars")
			, ifelse(aggregate.obj.hist, sum(nbins.hist > 0), sum(nbins.hist))
			, list(hist.ta = hist.ta.ref, hist.sl = hist.sl.ref, hist.tv = hist.tv.ref)
			, generations = generations, popsize = popsize
			, lower.bounds = attr(species.model, "lower.bounds")
			, upper.bounds = attr(species.model, "upper.bounds")
			, vectorized = FALSE, mprob = mprob
		)
	}
	
	if(!is.null(cl)) stopCluster(cl)
	attr(sol, "generations") <- generations
	attr(sol, "species.model") <- species.model
	return(sortSolutionParameters(sol))
}

## A convenience function to plot the evolution of the algorithm and assess convergence
generationPlot <- function(solutions, species.model = attr(solutions, "species.model")
	, plot.quantiles = c(0.10, 0.5, 0.90), only.pareto = FALSE
	, show.legend = TRUE, lwd = 1.5, mar = c(2.3, 2.3, 0.2, 2.3)
	, mgp = c(1.2, 0.2, 0), tcl = -0.25, ...) {

	generations <- attr(solutions, "generations")
	
	# some fancy colors and their semi-transparent counterparts
	plot.colors = matrix(
		c("#000000", "#00000022"
		, "#00bb00", "#00bb0022"
		, "#ff0000", "#ff000022"
		, "#ff9900", "#ff990022"
		, "#0077ff", "#0077ff22"
		, "#9900ff", "#99007722"), ncol = 2, byrow = T)
	plot.lty = rep(c(1, 2), each = nrow(plot.colors))
	plot.colors = rbind(plot.colors, plot.colors)	# duplicate
	
	# compute quantiles for plotting
	if(only.pareto) {
		quantiles <- vapply(solutions, function(x) {
			apply(x$par[x$pareto, , drop = FALSE], 2, quantile, plot.quantiles)
		}, matrix(0, nrow=3, ncol=attr(species.model, "npars")))
	} else {
		quantiles <- vapply(solutions, function(x) {
			apply(x$par[, , drop = FALSE], 2, quantile, plot.quantiles)
		}, matrix(0, nrow=3, ncol=attr(species.model, "npars")))
	}

	dimnames(quantiles) <- list(
		quantile = plot.quantiles
		, parameter = seq_len(attr(species.model, "npars"))
		, generation = generations)

	par(mar = mar, mgp = mgp, tcl = tcl, ...)
	plot.new()
	# plot scale for concentration and probabilites (left axis)
	plot.window(xlim = c(0, max(generations)), ylim = c(0, 1))
	axis(1)
	axis(2)
	# TODO: better way to tell apart the types of parameters
	coor.probs <- which(attr(species.model, "upper.bounds") <= 1)
	the.others <- setdiff(seq_len(attr(species.model, "npars")), coor.probs)
	for(p in coor.probs) {
# concentration and probabilities are those parameters bounded by 1
		polygon(c(generations, rev(generations)), c(quantiles[1, p, ]
			, rev(quantiles[3, p, ]))
			, border = NA, col = plot.colors[p, 2])
		lines(generations, quantiles[2, p, ], lwd = lwd, col = plot.colors[p, 1]
			, lty = plot.lty[p])
	}
	# plot scale for step lengths (right axis)
	plot.window(xlim = c(0, max(generations)), ylim = c(0, max(attr(species.model
		, "upper.bounds"))))
	axis(4)
	for(p in the.others) {		# step lengths are the last params
		polygon(c(generations, rev(generations)), c(quantiles[1, p, ]
			, rev(quantiles[3, p, ]))
			, border = NA, col = plot.colors[p, 2])
		lines(generations, quantiles[2, p, ], lwd = lwd, col = plot.colors[p, 1]
			, lty = plot.lty[p])
	}
	
	title(xlab = c("Generation"))
	mtext("Concentration and probabilities", 2, line = 1.2)
	mtext("Step lengths", 4, line = 1.2)
	box(bty = "u")

	if(show.legend) {
		legend("bottomleft"#, xpd = T, inset = c(-0, -0.15)
			, legend = attr(species.model, "param.names")
			, lwd = lwd, col = plot.colors[seq_len(attr(species.model, "npars")), 1]
			, lty = plot.lty[seq_len(attr(species.model, "npars"))]
			, box.lwd = 0, bg = "#ffffff77")
	}
	return(invisible(quantiles))
}

# Computes the turning angle variation histogram (variation computed in a
# time-moving window)
# This histogram gives a view of the relative amount of time that the animal
# spends in each distinct state
computeVariationHistogram <- function(relocs, nbins = 7, range = NULL, window.size = dim(relocs$stats)[1] %/% 100) {
	a.var <- angle.variation(relocs, window.size = window.size)
	if(is.null(range)) 
		range <- range(a.var)
	hist.var.ref <- binCounts(a.var, range, nbins)
	return(hist.var.ref)
}

angle.variation <- function(relocs, nlags = 100
	, window.size = dim(relocs$stats)[1] %/% nlags) {

	ta <- relocs$stats[, 2]
	sda <- numeric(length(ta) %/% window.size)

	for(i in seq_along(sda)) {
		sda[i] <- sd(ta[((i - 1) * window.size + 1) : (i * window.size)])
	}
	return(sda)
}

binCounts <- function(data, range, nbins, log = FALSE) {
	if(nbins == 0) return(NULL)
	if(log) {
		data = log(data + 1)
		range = log(range + 1)
	}
	bins <- seq(range[1], range[2], len = nbins + 1)
	inter <- findInterval(data, bins, rightmost.closed = TRUE)
	inter <- inter[inter > 0 & inter <= nbins]
	tr <- table(inter)
	tra <- rep(0, nbins)
	names(tra) <- paste(round(bins[seq_len(nbins)]), "-", round(bins[2:(nbins + 1)]), sep="")
	tra[as.numeric(names(tr))] <- tr
	return(tra)
}

# This function sorts the parameters of a population of solutions such that
# the states will be arranged with a decreasing order turning angle concentration
sortSolutionParameters <- function(solutions) {
	attrs <- attributes(solutions)
	spmodel <- attr(solutions, "species.model")
	if(inherits(solutions, "nsga2")) {
		return(sortSolutionParametersSingleGeneration(solutions, spmodel))
	} else if(inherits(solutions, "nsga2.collection")) {
		out <- lapply(solutions, sortSolutionParametersSingleGeneration, spmodel)
		attributes(out) <- attrs
		return(out)
	} else
		stop("Expecting an object of class nsga2")
}

sortSolutionParametersSingleGeneration <- function(solutions, spmodel) {
	if(!inherits(solutions, "nsga2")) stop("Expecting an object of class nsga2")
	
	types <- attr(spmodel, "param.types")
	npars <- attr(spmodel, "npars")
	correls <- grep("^TA", types)
	if(length(correls) > 9) stop("A maximum of 9 states is supported")

	types <- vapply(types, function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=""), "")

	tmp <- t(apply(solutions$par, 1, function(sol) {
		map <- order(sol[correls], decreasing=TRUE)
		newtypes <- types
		newtypes <- chartr(paste(seq_along(map), collapse=""), paste(map, collapse=""), newtypes)
		newtypes <- vapply(newtypes, function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=""), "")
		neworder <- match(newtypes, types)
		return(sol[neworder])
	}))
	solutions$par <- tmp
	return(solutions)
}

