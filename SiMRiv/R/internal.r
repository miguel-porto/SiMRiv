setClass("perceptionWindow",slots=c("type"="character","parameters"="numeric"))
setClass("state",slots=c("turningAngleConcentration"="numeric","perceptionWindow"="perceptionWindow","stepLength"="numeric","name"="character"))
setClass("species",slots=c("name"="character","transitionMatrix"="matrix","states"="list"))
#setClass("individual",slots=c("species"="species"))

MULTIPLIER=1000000L

perceptionWindow<-function(type = "circular", radius) {
	switch(pmatch(type,c("circular","gaussian"),nomatch=3),{
			out=new("perceptionWindow",type="circular",parameters=c(radius=radius))
		},{
			out=new("perceptionWindow",type="gaussian",parameters=c(sigma=radius))
		}
		,stop("'type' must be one of 'circular', 'gaussian'.")
	)
	return(out)
}

.printPercWind <- function(object) {
	switch(pmatch(object@type,c("circular","gaussian"),nomatch=3),{
		cat("Circular perception window with radius =",object@parameters[1],"\n")
	},{
		cat("Gaussian perception window with sigma =",object@parameters[1],"\n")
	})	
}

setMethod("show", signature(object="perceptionWindow"),.printPercWind)
setMethod("print", signature(x="perceptionWindow"),function(x) show(x))

# steplen is assumed to be in the same units of the coordinates! Usually meters.
state<-function(correlation, pwind = perceptionWindow("circular", 0), steplen = 1, name = "") {
	if(correlation==1) correlation=0.99999
	if(correlation==0) correlation=0.00001
	
	return(new("state"
		, turningAngleConcentration=correlation
		, perceptionWindow=pwind
		, stepLength=steplen
		, name=name ))
}

setValidity("state",function(object) {
	if(object@turningAngleConcentration<0 || object@turningAngleConcentration>1) return("Turning angle concentration must be within [0,1]")
	return(TRUE)
})


state.Resting<-function() {
	return(state(0,steplen=0,name="Resting"))
}

state.RW<-function() {
	return(state(0,steplen=1,name="Random Walk"))
}

state.CRW<-function(correlation) {
	return(state(correlation,steplen = 1,name = "Correlated Random Walk"))
}

.printState<-function(object) {
	if(object@name == "") {
		cat("Movement state with\n\t")
	} else {
		cat("State",object@name,"with\n\t")
	}
	show(object@perceptionWindow)
	cat("\tStep length =",object@stepLength,"\n\tTurning angle concentration =",object@turningAngleConcentration,"\n")
}

setMethod("show", signature(object="state"),.printState)

species<-function(states, trans = transitionMatrix(), name = "<unnamed>") {
	if(!(inherits(states,"list"))) states=list(states)

	out=new("species", states=states
		, transitionMatrix=trans
		, name=name
	)
	
	names=character(length(states))
	for(i in 1:length(names)) {
		if(states[[i]]@name == "") {
			names[i]=paste("State",i)
		} else {
			names[i]=paste(i,": ",states[[i]]@name,sep="")
		}
	}
	
	colnames(out@transitionMatrix)=names
	rownames(out@transitionMatrix)=names
	return(out)
}

.printSpecies<-function(object) {
	cat("Species",ifelse(length(object@name)==0,"<unnamed>",object@name),"with",length(object@states),"states:\n")
	for(i in 1:length(object@states)) {
		cat("\n",i,": ",sep="")
		show(object@states[[i]])
	}
	cat("\nState transition matrix\n")
	print(object@transitionMatrix)
}

setMethod("show", signature(object="species"), .printSpecies)

setValidity("species",function(object) {
	trans=object@transitionMatrix
	states=object@states
	if(!(inherits(states,"list"))) return("States must be a list of states")
	for(i in 1:length(states)) {
		if(!inherits(states[[i]],"state")) return("States must be a list of states")
	}
	if(dim(trans)[1]!=dim(trans)[2]) return("Transition matrix must be square")
	if(length(states)!=dim(trans)[1]) return("Transition matrix must have the same number of rows as the number of states provided")
	if(any(trans<0) || any(trans>1)) return("Transition matrix values must be probabilities between 0 and 1")
	if(any(apply(trans,1,sum)!=1)) return("Transition matrix must sum to 1 in all rows (not columns)")
	return(TRUE)
})

transitionMatrix <- function(...) {
	probs=unlist(list(...))
	if(any(probs<0) || any(probs>1)) stop("Probabilities must be between 0 and 1")
	N=(1+sqrt(1+4*length(probs)))/2
	if(round(N) != N) stop("You must specify all non-diagonal elements of the square matrix (by row)")
	out=matrix(nrow = N, ncol = N)
	c=1
	for(i in 1:N) {
		for(j in 1:N) {
			if(i == j) next
			out[i,j] = probs[c]
			c = c + 1
		}
	}
	diag(out) = 1 - apply(out,1,sum,na.rm=T)
	if(any(apply(out,1,sum) != 1)) return("Transition matrix must sum to 1 in all rows (not columns)")
	return(out)
}
#individual<-function(species) {
#	return(new("individual",species=species))
#}

#.printIndividual<-function(object) {
#	cat("Individual of species",object@species@name,"\n")
#}

#setMethod("show", signature(object="individual"),.printIndividual)

