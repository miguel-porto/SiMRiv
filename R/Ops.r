setClass("perceptualRange",slots=c("type"="character","parameters"="numeric"))
setClass("state",slots=c("turningAngleConcentration"="numeric","perceptualRange"="perceptualRange","stepLength"="numeric","name"="character"))
setClass("species",slots=c("name"="character","transitionMatrix"="matrix","states"="list", "resistanceMap" = "ANY"))
#setClass("individual",slots=c("species"="species"))

setMethod("+", signature(e1 = "species"), function(e1, e2) {
	for(i in 1:length(e1@states)) {
		e1@states[[i]] <- callGeneric(e1@states[[i]], e2)
	}
	return(e1)
})

setMethod("*", signature(e1 = "species"), function(e1, e2) {
	for(i in 1:length(e1@states)) {
		e1@states[[i]] <- callGeneric(e1@states[[i]], e2)
	}
	return(e1)
})

setMethod("+", signature(e1 = "state", e2 = "numeric"), function(e1, e2) {
	e1@stepLength <- e2
	return(e1)
})

setMethod("*", signature(e1 = "state", e2 = "numeric"), function(e1, e2) {
	e1@perceptualRange@parameters[["radius"]] = e2
	return(e1)
})

# Combine states for multistate movements
setMethod("+", signature(e1 = "state", e2 = "state"), function(e1, e2) {
	return(list(e1, e2))
})

setMethod("+", signature(e1 = "list", e2 = "state"), function(e1, e2) {
	return(c(e1, e2))
})

setMethod("+", signature(e1 = "state", e2 = "list"), function(e1, e2) {
	return(c(e1, e2))
})

