simulate <- function(individuals, time, coords = NULL, states = NULL, resist = NULL) {
	if(mode(time)!="numeric") stop("time must be numeric")
	if(!is.null(coords)) {
		if(!inherits(coords,"matrix") || dim(coords)[2]!=2) {
			if(!inherits(coords, "numeric") || length(coords) != 2) {stop("coords must be a 2-column matrix with initial coordinates or a vector of length 2")}
			coords=cbind(as.integer(coords[1]),as.integer(coords[2]))
		} else {coords=cbind(as.integer(coords[,1]),as.integer(coords[,2]))}
	} else {
		if(missing(resist) || is.null(resist)) {
			coords=cbind(
				as.integer(rep(0,length(individuals)))
				,as.integer(rep(0,length(individuals)))
			)
			warning("No starting coordinates and no raster given: starting positions set to (0,0)")
		} else {
			e=extent(resist)
			coords=cbind(
				as.integer(runif(length(individuals),e@xmin,e@xmax))
				,as.integer(runif(length(individuals),e@ymin,e@ymax))
			)
		}
	}

	if(!is.null(states)) {
		stop("Initial states not implmented yet")
		if(mode(states)!="numeric" || length(states)!=length(individuals)) stop("states must be an integer vector of initial states, the same length of individuals")
		# TODO: initial states not implmented
	}
	
	if(!(inherits(individuals,"list"))) {
		if(!(inherits(individuals,"species")))
			stop("individuals must be a list of species class")
		else
			individuals=list(individuals)
	}

	if(dim(coords)[1] != length(individuals)) stop("The number of rows in the 'coords' matrix must be the same as the number of individuals")
		
	for(i in 1:length(individuals)) {
		if(!inherits(individuals[[i]],"species")) stop("individuals must be a list of species class")
	}
	
	.Call("_simulate_individuals",individuals,coords,as.integer(time),resist,new.env())
}


#resistanceFromSHP<-function(shp,res,value=1) {
#	s=readShapeSpatial(shp)
#	a=raster(ext=extent(s),crs=proj4string(s),resolution=res)
#	b=rasterize(s,a,value)
#	return(b)
#}

sampleMovement<-function(relocs, resolution = 1, resist = NULL) {
	if(resolution<1) stop("Resolution must be at least 1 time tick.")
	tmp=as.integer(round(resolution))
	if(tmp != resolution) stop("Resolution must be an integer number")
	if(dim(relocs)[2] != 3) stop("Only implemented for single individual simulations at the moment.")
	relocs=relocs[seq(1,dim(relocs)[1],by=tmp),1:2]
	diffs=apply(relocs,2,diff)
	steplengths=sqrt(apply(diffs^2,1,sum))
	angles=atan2(diffs[,2],diffs[,1])
	turnangles=diff(angles)
	turnangles[turnangles>pi]=-2*pi+turnangles[turnangles>pi]
	turnangles[turnangles<(-pi)]=2*pi+turnangles[turnangles<(-pi)]
	#hist(turnangles)
	#plot(density(turnangles))
	if(!is.null(resist)) {
		resistance=.Call("stepRasterAccumulator",relocs,resist,new.env())
		stats=data.frame(
			steplengths=steplengths[2:length(steplengths)]
			,turningangles=turnangles
			,resistance=resistance[1:(length(resistance)-1)])
	} else {
		stats=data.frame(
			steplengths=steplengths[2:length(steplengths)]
			,turningangles=turnangles)
	}
	return(list(
		relocs=relocs
		,stats=stats
	))
}

