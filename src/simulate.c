// TODO: think better about how to use OMP better, we don't use it for now.
//#define USEOPENMP
#ifdef USEOPENMP
#include <omp.h>
#endif
#include "SiMRiv.h"
SEXP rho;

/**
* Main simulation function
*/
SEXP _simulate_individuals(SEXP _individuals, SEXP _starting_positions
	, SEXP _timespan, SEXP _angles, SEXP _resist, SEXP envir) {
	
	#ifdef USEOPENMP
	omp_set_num_threads(omp_get_num_procs ( ));
	#endif
	rho=envir;
	RASTER *resist = NULL;
	GetRNGstate();
	
	if(_resist != R_NilValue) {
		resist = openRaster(_resist,rho);
	}
	
	int timespan = INTEGER_POINTER(_timespan)[0]
		, *start = INTEGER_POINTER(_starting_positions);
	double *angles = NULL;
	double *prelocs;
	unsigned int ninds = LENGTH(_individuals), i, j, k, time, tmp1, tmp2;
	SEXP relocs, tmp3, tmp4;
	const char *tmp5;
	float curangtrans;
	
	if(_angles != R_NilValue)
		angles = NUMERIC_POINTER(_angles);

// pointers to individual data
	INDIVIDUAL *ind = malloc(sizeof(INDIVIDUAL) * ninds);
	
// output is a matrix with columns x,y,state appended for each individual	
	PROTECT(relocs = allocMatrix(REALSXP, timespan, 3 * ninds));
	prelocs=NUMERIC_POINTER(relocs);
	
// get pointers to individual and state data
	for(i = 0; i < ninds; i++) {
		ind[i].pspecies = VECTOR_ELT(_individuals, i);
		ind[i].transitionmatrix = NUMERIC_POINTER(GET_SLOT(ind[i].pspecies
			, SCALARCHAR("transitionMatrix")));
		ind[i].nstates = LENGTH(GET_SLOT(ind[i].pspecies, SCALARCHAR("states")));
		ind[i].states = malloc(sizeof(STATE) * ind[i].nstates);		
// pointers for each state data
		for(j = 0; j < ind[i].nstates; j++) {
			tmp3 = VECTOR_ELT(GET_SLOT(ind[i].pspecies, SCALARCHAR("states")), j);
			ind[i].states[j].taconc = NUMERIC_POINTER(GET_SLOT(tmp3
				, SCALARCHAR("turningAngleConcentration")))[0];
			ind[i].states[j].steplength = NUMERIC_POINTER(GET_SLOT(tmp3
				, SCALARCHAR("stepLength")))[0];
			tmp4 = GET_SLOT(tmp3, SCALARCHAR("perceptualRange"));
			ind[i].states[j].pwind.radius = NUMERIC_POINTER(GET_SLOT(tmp4
				, SCALARCHAR("parameters")))[0];
			tmp5 = CHAR(STRING_ELT(GET_SLOT(tmp4, SCALARCHAR("type")), 0));
			if(strcmp(tmp5, "circular") == 0)
				ind[i].states[j].pwind.type = CIRCULAR;
			else if(strcmp(tmp5,"gaussian") == 0)
				ind[i].states[j].pwind.type = GAUSSIAN;
			else error("Invalid perceptual range type.");
			
// let's create a lookup table for perception window weights!
			if(ind[i].states[j].pwind.type == GAUSSIAN) {
				float pwsigma = ind[i].states[j].pwind.radius;
				// NOTE: radius is taken as sigma and increased 4 times
				ind[i].states[j].pwind.radius =
					(int) ceil(ind[i].states[j].pwind.radius * 4);
				ind[i].states[j].pwind.weights =
					calloc(ind[i].states[j].pwind.radius + 1, sizeof(float));
				for(k = 0; k < ind[i].states[j].pwind.radius; k++) {
					ind[i].states[j].pwind.weights[k] =
						exp(-(double) k * k / (2 * pwsigma * pwsigma));
				}
			}
		}
	}

// START SIMULATION
	{
		PDF tmpPDF, tmprotPDF;
		CDF tmpMultCDF;
		double **curtrans = malloc(sizeof(double*) * ninds);
		unsigned long r, s, k;
		float lengthmove;
		STATE *tmpstate;
// assign initial states, angles and starting positions
		for(i = 0; i < ninds; i++) {
			ind[i].curpos.x = start[i];
			ind[i].curpos.y = start[i + ninds];
// random initial state
			ind[i].curstate = runif(0, ind[i].nstates - 1);
// uniform random angle if angle not provided
			ind[i].curang = angles ? (ISNAN(angles[i]) ? drawRandomAngle(NULL)
				: ((angles[i] + PI) / ANGLESTEP)) : drawRandomAngle(NULL);
// initialize states
			for(j = 0; j < ind[i].nstates; j++) {
// compute base circular PDFs for all states of all individuals (centered on 0),
// with given concentration
				circNormal(ind[i].states[j].taconc, ind[i].states[j].basePDF
					, ind[i].states[j].scaledPDF);
// compute the respective cumulative circular PDFs
				for(k = 1, ind[i].states[j].cumPDF[0] =
					(long)(ind[i].states[j].basePDF[0] * MULTIPLIER)
					; k < ANGLERES; k++) {
					ind[i].states[j].cumPDF[k] =
						ind[i].states[j].cumPDF[k - 1]
						+ (long) (ind[i].states[j].basePDF[k] * MULTIPLIER);
				}
			}
		}
/***************************************
** MAIN TIME LOOP
***************************************/
// NOTE: time is unitless for now
		for(time = 0; time < timespan; time++) {
			for(i = 0; i < ninds; i++)
// for now, constant transition matrix
				curtrans[i] = ind[i].transitionmatrix;

// LOOP FOR EACH INDIVIDUAL
// TODO: we can parallelize this if we just want a batch of non-interacting simulations
// tmp2 is just a relative pointer to output matrix
			for(i = 0, tmp2 = 0; i < ninds; i++, tmp2 += timespan) {
// draw new state according to transition matrix
				r = runif(0, MULTIPLIER - 1);
				for(k = 0, s = (unsigned long) (curtrans[i][ind[i].curstate]
						* MULTIPLIER)
					; k < (ind[i].nstates - 1) && r >= s
					; k++, s += (unsigned long) (curtrans[i][ind[i].curstate 
						+ k * ind[i].nstates] * MULTIPLIER)) {}

				ind[i].curstate = k;
// just a handy pointer to state data				
				tmpstate = &ind[i].states[ind[i].curstate];
// rotate the 0-centered base PDF (the one calculated from the state's
// concentration parameter) to be centered on the previous step angle
				rotatePDF(tmpstate->scaledPDF, tmprotPDF
					, ind[i].curang - ANGLECENTER);
// compute the empirical resistance PDF around the current position
				computeEmpiricalResistancePDF(ind[i].curpos, resist
					, &tmpstate->pwind, tmpPDF);
				if(tmpPDF[0] != -1) {		// resistance is heterogeneous
// multiply rotated base PDF by empirical PDF and compute cumulative PDF on the fly
					for(j = 1, tmpMultCDF[0] = (unsigned long) (tmprotPDF[0]
						* tmpPDF[0] * MULTIPLIER); j < ANGLERES; j++) {
						tmpMultCDF[j] = tmpMultCDF[j - 1] + (long) (tmprotPDF[j]
							* tmpPDF[j] * MULTIPLIER);
					}

// TODO: what to do when the desired direction is facing towards an
// infinite resistance area and there is no overlap of PDFs?
// Here we just draw a uniform random angle
					if(tmpMultCDF[ANGLERES - 1] == 0)
						ind[i].curang = drawRandomAngle(NULL);
					else
// otherwise, follow normal procedure: draw random angle based on compound PDF
// (resistance x correlated component)
						ind[i].curang = drawRandomAngle(tmpMultCDF);
				} else {	// resistance is equal in all directions, ignore
					for(j = 1, tmpMultCDF[0] = (unsigned long) (tmprotPDF[0]
						* MULTIPLIER); j < ANGLERES; j++)
						tmpMultCDF[j] = tmpMultCDF[j - 1] + (long) (tmprotPDF[j]
							* MULTIPLIER);
// draw random angle based only on base PDF (correlated component)
					ind[i].curang = drawRandomAngle(tmpMultCDF);
				}

// compute the realized length of this step, taking into account resistance
				curangtrans = ind[i].curang * ANGLESTEP - PI;
				lengthmove = computeLengthMove(tmpstate->steplength
					, ind[i].curpos, resist, curangtrans);
// take the move
				if(lengthmove > 0) {
					ind[i].curpos.x += cos(curangtrans) * lengthmove;
					ind[i].curpos.y += sin(curangtrans) * lengthmove;
				}
				
				if(isnan(ind[i].curpos.x)) {	// this should never happen, but...
					Rprintf("%f %f %f %f %f %f\n", tmpstate->steplength
						, ind[i].curpos.x, ind[i].curpos.y, lengthmove
						, ind[i].curang, curangtrans);
					Rf_error("Unexpected error, please report.");
				}

// write position and state in output matrix
// this is the column offset of the current individual in the output matrix
				tmp1 = tmp2 * 3;
				prelocs[tmp1 + time] = ind[i].curpos.x;
				prelocs[tmp1 + timespan + time] = ind[i].curpos.y;
				prelocs[tmp1 + timespan + timespan + time] = ind[i].curstate;
			}
		}
		
		free(curtrans);
	}

// free things	
	PutRNGstate();
	UNPROTECT(1);

	for(i = 0; i < ninds; i++) {
		for(j = 0; j < ind[i].nstates; j++) {
			if(ind[i].states[j].pwind.type == GAUSSIAN)
				free(ind[i].states[j].pwind.weights);
		}
		free(ind[i].states);
	}
	free(ind);
	if(_resist != R_NilValue) closeRaster(resist);
	return relocs;
}

/**
* Make a circular wrapped normal PDF centered on zero
* Note that there's an offset of -PI, so the first value corresponds to the
* density of -PI. If scaledout is provided, also output range-standardized PDF
* This code was adapted from the CircStats package, function dwrpnorm()
* _rho is the concentration parameter, varying from 0 (flat) to 1 (one peak)
*/
void circNormal(float _rho, float* out, float* scaledout) {
	float var = -2.f * log(_rho), next, last, delta;
	float ang, max = -1;
	int k, i;
    
	for(i = 0, ang = -PI; i <= ANGLECENTER; i++, ang += ANGLESTEP) {
		k = 0;
		next = circNormalTerm(ang, var, k);
    	delta = 1;
		while(delta > TOLERANCECIRCNORM) {
		    k++;
		    last = next;
		    next = last + circNormalTerm(ang, var, k)
		    	+ circNormalTerm(ang, var, -k);
		    delta = fabs(next - last);
		}
		if(next > max) max = next;
		out[i] = next;
	}
	for(i = ANGLECENTER; i < ANGLERES; i++) out[i] = out[(ANGLERES - 1) - i];
	
	if(scaledout) {		// range-standardize
		for(i=0; i < ANGLERES; i++) scaledout[i] = out[i] / max;
	}
}

/**
* Draws a random angle from CDF or a uniform random angle if cdf==NULL
*/
float drawRandomAngle(CDF cdf) {
	if(cdf == NULL) {
		return runif(0, ANGLERES - 1);
	} else {
		if(cdf[ANGLERES-1] == 0) return runif(0, ANGLERES - 1);
		return densityRand(ANGLERES, cdf);
	}
}

/**
* Draws a random number from a discrete CDF
*/
int densityRand(int nValues, unsigned long *cdf) {
	long lo, hi, mid;
	float r;
	lo = 0;
	hi = nValues - 1;
	r = cdf[hi] * unif_rand();

// binary search
	while(lo < hi) {
		mid = (lo + hi) / 2;
		if (r < cdf[mid]) hi = mid; else lo = mid + 1;
	}
	return lo;
}

/**
* Computes the actual distance moved in a given direction, taking into account
* the resistance that has to be crossed between the two desired points.
* NOTE: currently, this only takes into account the resistance of the starting
* and projected ending positions (average both values). But we can improve this,
* by taking into account also some intermediate positions. To what extent would
* such improvement degrade performance?
*/
inline float computeLengthMove(double baseStepLength, POINT curpos
	, RASTER *resist, float angle) {
	if(resist == NULL) return baseStepLength;
// here we take a simple approach: just compute the mean resistance between
// the starting point and the (theoretical) end point
	float tmp1 = (float) extractRasterValueNoNaN(resist, curpos.x + cos(angle)
		* baseStepLength, curpos.y + sin(angle) * baseStepLength);
// if it is infinite resistence, do not go there at all!
// NOTE: should we leave this condition here or just let it go?
	if(tmp1 > 0.999) return 0;
	float tmp = (float) extractRasterValueNoNaN(resist, curpos.x, curpos.y);
// step length will be proportional to the mean "conductance"
	return baseStepLength * (1 - (tmp + tmp1) / 2);
}

/*
* Computes the circular EPDF at a given location, according to resistance values
* in the surroundings (along discrete radial lines), and angular bias.
* Returns the PDF in the provided pointer.
* NOTE: the PDF has an arbitrary scale, the integral is not constant!! We don't
* need it to be, here.
*/
void computeEmpiricalResistancePDF(POINT curpos,
	const RASTER *resist,PERCEPTIONWINDOW *percwind,PDF pdf) {

	bool allinf = true;
	float step;
	
	if(resist == NULL) {	// signal that a uniform random should be drawn
		pdf[0] = -1;
		return;
	}
	
	switch(percwind->type) {
	case CIRCULAR:
		step = percwind->radius / ACCUMULATORRESOLUTION;
		int i, j;
		float tcos, tsin, ang, sum;
		POINT tmppos;
		double tmp;
		
		#ifdef USEOPENMP
		#pragma omp parallel for firstprivate(step,resist,curpos,percwind) private(i,j,tcos,tsin,ang,sum,tmppos,tmp) shared(allinf,pdf)
		#endif
		for(i = 0; i < ANGLERES; i++) {	// make a whole circle in radial lines
			ang = -PI + i * ANGLESTEP;
			tmppos = curpos;
// for each angle, sum the resistance values along the radial line centered on
// curpos
// TODO: create a look-up table to speed up things here? Since ang is discrete
			tcos = cos(ang) * step;
			tsin = sin(ang) * step;
// integrate conductance along each radial line
			for(j = 0, sum = 0; j < percwind->radius; j += step) {
				tmp = extractRasterValue(resist, tmppos.x, tmppos.y);

				if(isnan(tmp) || tmp == NA_REAL || tmp == 1)
					sum += 0;
				else {
					sum += 1 - tmp;
					allinf = false;
				}
				tmppos.x += tcos;
				tmppos.y += tsin;
			}
			pdf[i] = sum;
		}
		break;
		
	case GAUSSIAN:	// gaussian perceptual range (i.e. with gaussian weights)
// TODO is this adequate in the gaussian?
		step = percwind->radius / ACCUMULATORRESOLUTION;
		#ifdef USEOPENMP
		#pragma omp parallel
		#endif
		{
			int i, j;
			float tcos, tsin, ang, sum;
			POINT tmppos;
			double tmp;
			
			#ifdef USEOPENMP
			#pragma omp for
			#endif
			for(i=0; i < ANGLERES; i++) {	// make a whole circle
				ang = -PI + i * ANGLESTEP;
				tmppos = curpos;
				tcos = cos(ang) * step;
				tsin = sin(ang) * step;
// for each angle, sum the resistance values along the radial line centered on
// curpos
				for(j = 0, sum = 0; j < percwind->radius; j += step) {
					tmp = extractRasterValue(resist, tmppos.x, tmppos.y);

					if(isnan(tmp) || tmp == NA_REAL || tmp == 1)
						sum += 0;
					else {
						sum += (1 - tmp) * percwind->weights[(int) floor(j)];
						allinf = false;
					}
					tmppos.x += tcos;
					tmppos.y += tsin;
				}
				
				pdf[i] = sum;
			}
		}
		break;
		
	default:
		Rprintf("PW: %d\n", (int)percwind->type);
		error("Perceptual range type not implemented yet");
		break;
	}
	if(allinf) pdf[0] = -1;
}

/**
* Shifts the circular PDF by and angle offset
*/
void rotatePDF(PDF pdf, PDF out, float ang) {
	int offset = ANGLERES - ang, j, i;
	for(i = 0, j = offset; i < ANGLERES; i++, j++)
		out[i] = pdf[(j < 0 ? (ANGLERES + j) : (j % ANGLERES))];
}

/**
* Computes the mean value of the given raster "infinitesimally" along each step
* of the movement. I.e., integrates the raster along each segment of the path
* and divides by its length.
* relocs is assumed to be a matrix whose first 2 columns are X and Y.
* "infinitesinally" means that it discretizes each step into
* ACCUMULATORRESOLUTION pieces, so the microstep length is variable.
*/
SEXP stepRasterAccumulator(SEXP relocs, SEXP _resist, SEXP envir) {
	if(_resist == R_NilValue) error("Must provide a raster");
	rho = envir;
	RASTER *resist;
	SEXP out;
	double *prelocs = NUMERIC_POINTER(relocs), tmp, *pout;
// -1 cause the loop excludes the last one
	int nrelocs = (int) INTEGER_POINTER(GET_DIM(relocs))[0] - 1;
	double vx, vy;
	float sum;
	POINT tmppos, endpoint;
	int i, j, nsteps;
	
	resist = openRaster(_resist, rho);
// output is a vector with the accumulated resistance in each step
	PROTECT(out = NEW_NUMERIC(nrelocs));
	pout = NUMERIC_POINTER(out);
// loop of the steps
	for(i = 0; i < nrelocs; i++) {
		tmppos.x = prelocs[i];
		tmppos.y = prelocs[i + nrelocs + 1];
		endpoint.x = prelocs[i + 1];
		endpoint.y = prelocs[i + nrelocs + 2];
// compute a vector with length of microstep, with the direction of the step
		vx = endpoint.x - tmppos.x;
		vy = endpoint.y - tmppos.y;
		vx /= ACCUMULATORRESOLUTION;
		vy /= ACCUMULATORRESOLUTION;
		nsteps = ACCUMULATORRESOLUTION;
		
// "integrate" raster values between starting and ending points of this step
		for(j = 0, sum = 0; j < nsteps; j++, tmppos.x += vx, tmppos.y += vy) {
// TODO: we could optimize this, we don't need to call a function with calculations
// here, just integrate and avoid calculations
			tmp = extractRasterValue(resist, tmppos.x, tmppos.y);
			if(isnan(tmp) || tmp == NA_REAL || tmp == 1)
				sum += 1;
			else
				sum += tmp;
		}
		pout[i] = sum / nsteps;
	}
	
	closeRaster(resist);
	UNPROTECT(1);
	return out;
}

