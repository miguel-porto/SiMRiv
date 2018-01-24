MULTIPLIER <- 1000000L

perceptualRange <- function(type = "circular", radius) {
	switch(pmatch(type, c("circular","gaussian"), nomatch=3), {
			out <- new("perceptualRange", type="circular"
				, parameters = c(radius=radius))
		},{
			out <- new("perceptualRange", type="gaussian"
				, parameters=c(sigma=radius))
		}
		, stop("'type' must be one of 'circular', 'gaussian'.")
	)
	return(out)
}

.printPercWind <- function(object) {
	switch(pmatch(object@type, c("circular","gaussian"), nomatch=3), {
		cat("Circular perceptual range with radius =", object@parameters[1], "\n")
	}, {
		cat("Gaussian perceptual range with sigma =", object@parameters[1], "\n")
	})	
}

setMethod("show", signature(object="perceptualRange"), .printPercWind)
setMethod("print", signature(x="perceptualRange"), function(x) show(x))

# steplen is assumed to be in the same units of the coordinates! Usually meters.
state <- function(concentration, pwind = perceptualRange("circular", 0)
	, steplen = 1, name = "") {

	if(concentration == 1) concentration <- 0.99999
	if(concentration == 0) concentration <- 0.00001
	
	return(new("state"
		, turningAngleConcentration=concentration
		, perceptualRange=pwind
		, stepLength=steplen
		, name=name ))
}

setValidity("state", function(object) {
	if(object@turningAngleConcentration < 0
		|| object@turningAngleConcentration > 1)
		return("Turning angle concentration must be within [0,1]")
	return(TRUE)
})
 
state.Resting <- function() {
	return(state(0, steplen = 0, name = "Rest"))
}

state.RW <- function() {
	return(state(0, steplen = 1, name = "RW"))
}

state.CRW <- function(concentration) {
	return(state(concentration, steplen = 1, name = "CRW"))
}

.printState <- function(object) {
	if(object@name == "") {
		cat("Movement state with\n\t")
	} else {
		cat("State",object@name,"with\n\t")
	}
	show(object@perceptualRange)
	cat("\tMax step length =", object@stepLength
		, "\n\tTurning angle concentration =", object@turningAngleConcentration
		, "\n")
}

setMethod("show", signature(object="state"), .printState)

species <- function(states, trans = transitionMatrix(), name = "<unnamed>"
	, resistanceMap = NULL) {
	if(!(inherits(states, "list"))) states <- list(states)

	out=new("species", states = states
		, transitionMatrix = trans
		, name = name
		, resistanceMap = resistanceMap
	)
	
	names=character(length(states))
	for(i in seq_along(names)) {
		if(states[[i]]@name == "") {
			names[i] <- paste("State", i)
		} else {
			names[i] <- paste(i, ": ", states[[i]]@name, sep="")
		}
	}
	
	colnames(out@transitionMatrix) <- names
	rownames(out@transitionMatrix) <- names
	return(out)
}

.printSpecies <- function(object) {
	cat("Species", ifelse(length(object@name) == 0, "<unnamed>", object@name)
		, "with", length(object@states), "states:\n")
	for(i in seq_along(object@states)) {
		cat("\n", i, ": ", sep="")
		show(object@states[[i]])
	}
	cat("\nState transition matrix\n")
	print(round(object@transitionMatrix, 4))
	#cat("\nResistance map\n")
	#print(object@resistanceMap)
}

setMethod("show", signature(object="species"), .printSpecies)

setValidity("species", function(object) {
	trans <- object@transitionMatrix
	states <- object@states
	if(!(inherits(states, "list"))) return("States must be a list of states")
	for(i in seq_along(states)) {
		if(!inherits(states[[i]], "state"))
			return("States must be a list of states")
	}
	if(dim(trans)[1] != dim(trans)[2]) return("Transition matrix must be square")
	if(length(states) != dim(trans)[1])
		return("Transition matrix must have the same number of rows as the number of states provided")
	if(any(trans<0) || any(trans>1))
		return("Transition matrix values must be probabilities between 0 and 1")
	if(any(abs(apply(trans, 1, sum) - 1) > 0.000001 ))
		return("Transition matrix must sum to 1 in all rows (not columns)")
	return(TRUE)
})

transitionMatrix <- function(...) {
	probs <- unlist(list(...))
	if(any(probs < 0) || any(probs > 1))
		stop("Probabilities must be between 0 and 1")
	N <- (1 + sqrt(1 + 4 * length(probs))) / 2
	if(round(N) != N)
		stop("You must specify all non-diagonal elements of the square matrix (by row)")
	out <- matrix(nrow = N, ncol = N)
	c <- 1
	for(i in seq_len(N)) {
		for(j in seq_len(N)) {
			if(i == j) next
			out[i, j] = probs[c]
			c = c + 1
		}
	}
	diag(out) = 1 - apply(out, 1, sum, na.rm=TRUE)
	#if(any(apply(out,1,sum) != 1)) stop("Transition matrix must sum to 1 in all rows (not columns)")
	return(out)
}

