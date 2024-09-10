simulate <- function(individuals, time, coords = NULL, states = NULL
	, resist = NULL, angles = NULL, start.resistance) {
# TODO: different resistance raster for each individual,
# using the species' resistanceMap
	if(mode(time) != "numeric") stop("time must be numeric")

	if(!(inherits(individuals, "list"))) {
		if(!(inherits(individuals, "species")))
			stop("individuals must be a list of species class")
		else
			individuals <- list(individuals)
	}

	if(!is.null(coords)) {
		if(!inherits(coords, "matrix") || dim(coords)[2] != 2) {
			if(!inherits(coords, "numeric") || length(coords) != 2) 
				stop("coords must be a 2-column matrix with initial coordinates or a vector of length 2")
			coords <- cbind(as.integer(coords[1]), as.integer(coords[2]))
		} else {
			coords <- cbind(as.integer(coords[,1]), as.integer(coords[, 2]))
		}
	} else {
		if(missing(resist) || is.null(resist)) {
			coords <- cbind(
				as.integer(rep(0, length(individuals)))
				, as.integer(rep(0, length(individuals)))
			)
# warning("No starting coordinates and no raster given: starting positions set to (0,0)")
		} else {
			if(missing(start.resistance) || is.null(start.resistance)) {
				e <- as.vector(terra::ext(resist))
				coords <- cbind(
					as.integer(runif(length(individuals), e["xmin"], e["xmax"]))
					, as.integer(runif(length(individuals), e["ymin"], e["ymax"]))
				)
			} else {
				coords <- matrix(as.integer(
					terra::xyFromCell(resist, sample(which(terra::values(resist) <= start.resistance)
					, length(individuals)))
				), ncol = 2)
			}
		}
	}

	if(!is.null(states)) {
		stop("Initial states not implmented yet")
		if(mode(states)!="numeric" || length(states) != length(individuals))
			stop("states must be an integer vector of initial states, the same length of individuals")
		# TODO: initial states not implmented
	}
	
	if(!is.null(angles)) {
		if(length(angles) != length(individuals))
			stop("angles must be a numeric vector of initial directions in radians, the same length of individuals")
		angles <- (-angles + pi / 2) %% (2 * pi)
	}

	if(dim(coords)[1] != length(individuals))
		stop("The number of rows in the 'coords' matrix must be the same as the number of individuals")
		
	for(i in seq_along(individuals)) {
		if(!inherits(individuals[[i]], "species"))
			stop("individuals must be a list of species class")
	}
	.Call(SR__simulate_individuals, individuals, coords, as.integer(time), angles, resist, new.env())
}

resistanceFromShape <- function(shp, baseRaster, res, binary = is.na(field)
	, field = NA, background = 1, buffer = NA, margin = 0, mapvalues = NA
	, extend = TRUE, ...) {

	if(missing(baseRaster) && missing(res))
		stop("Either raster resolution or a base raster must be given")
	if(inherits(shp, "character")) {
		l <- sf::st_read(shp)
	} else {
		l <- sf::st_as_sf(shp)
	}

	if(!all(is.na(buffer))) {
		b <- sf::st_buffer(l, dist = buffer)
	} else {
		b <- l
	}

	if(missing(baseRaster)) {
		er <- terra::rast(ext = terra::ext(b) + margin, crs = sf::st_crs(l)$proj4string
			, resolution = res)
	} else {
		if(extend) {
			er <- extend(baseRaster, terra::ext(b) + margin, fill = background)
		} else {
			er <- baseRaster
		}
	}
	
	if(binary) {
		r <- rasterize(b, er, background = 1, field = 0
			, update = !missing(baseRaster))
	} else {
		if(inherits(field, "numeric")) {
			r <- rasterize(b, er, field = field, background = background
				, update = !missing(baseRaster), ...)
		} else {
			if(all(is.na(mapvalues))) {
				if(!inherits(b[, field, drop = TRUE], "numeric"))
					stop("field must be a numeric field in the interval [0, 1], otherwise you have to specify mapvalues for translating field values")
				tmp <- as.numeric(b[, field, drop=TRUE])
				tmp[is.na(tmp)] <- background
				b[, field] <- tmp
				r <- rasterize(b, er, field = field, background = background
					, update = !missing(baseRaster), ...)
			} else {
				if(!inherits(mapvalues, "numeric"))
					stop("mapvalues must be a named numeric vector in the interval [0, 1]")
				empty <- names(mapvalues) == ""
				if(any(empty)) {
					emptyvalue <- mapvalues[empty]
					mapvalues <- mapvalues[!empty]
				}
				tmp <- mapvalues[as.character(b[, field, drop=TRUE])]
				names(tmp) <- NULL
				if(exists("emptyvalue")) {
					tmp[is.na(tmp)] <- emptyvalue
				} else {
					tmp[is.na(tmp)] <- background
				}
				b[, field] <- tmp
				r <- rasterize(b, er, field = field, background = background
					, update = !missing(baseRaster), ...)
			}
		}
	}

	if(terra::global(r, min) < 0 || terra::global(r, max) > 1)
		warning("Resistance values must be in the interval [0, 1]. Use 'mapvalues' to fix this.")
		
	return(r)
}

sampleMovement <- function(relocs, resolution = 1, resist = NULL) {
	if(resolution < 1) stop("Resolution must be at least 1 time tick.")
	tmp <- as.integer(round(resolution))
	if(tmp != resolution) stop("Resolution must be an integer number")
	if(dim(relocs)[2] > 3)
		warning("Only implemented for single individual simulations at the moment, using only the first two columns as coordinates")
	relocs <- relocs[seq(1, dim(relocs)[1], by=tmp), 1:2]
	diffs <- apply(relocs, 2, diff)
	steplengths <- sqrt(apply(diffs ^ 2, 1, sum))
	angles <- atan2(diffs[, 2], diffs[, 1])
	turnangles <- diff(angles)
	nas <- is.na(turnangles)
	turnangles[nas] <- 0
	turnangles[turnangles > pi] <- -2 * pi + turnangles[turnangles > pi]
	turnangles[turnangles < -pi] <- 2 * pi + turnangles[turnangles < -pi]
	turnangles[nas] <- NA

	if(!is.null(resist)) {
		resistance <- .Call(SR_stepRasterAccumulator, relocs, resist, new.env())
		stats <- data.frame(
			steplengths = steplengths[2:length(steplengths)]
			, turningangles = turnangles
			, resistance = resistance[1:(length(resistance) - 1)])
	} else {
		stats <- data.frame(
			steplengths = steplengths[2:length(steplengths)]
			, turningangles = turnangles)
	}
	return(list(
		relocs = relocs
		, stats = stats
	))
}

