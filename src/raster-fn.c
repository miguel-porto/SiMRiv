#include "SiMRiv.h"

/**
* A couple of functions imported from the raster package
*/
SEXP getRasterExtent(SEXP raster,SEXP rho) {
	SEXP ans,s,t;
	if(!isEnvironment(rho)) error("'rho' should be an environment");

	t = s = PROTECT(allocList(2));
	SET_TYPEOF(s, LANGSXP);
	SETCAR(t, install("extent")); t = CDR(t);
	SETCAR(t, raster);
	ans=eval(s, rho);
	UNPROTECT(1);
	return ans;
}

SEXP getRasterDim(SEXP raster,SEXP rho) {
	SEXP ans,s,t;
	if(!isEnvironment(rho)) error("'rho' should be an environment");

	t = s = PROTECT(allocList(2));
	SET_TYPEOF(s, LANGSXP);
	SETCAR(t, install("dim")); t = CDR(t);
	SETCAR(t, raster);
	ans=eval(s, rho);
	UNPROTECT(1);
	return ans;
}

SEXP getRasterRes(SEXP raster,SEXP rho) {
	SEXP ans,s,t;
	if(!isEnvironment(rho)) error("'rho' should be an environment");

	t = s = PROTECT(allocList(2));
	SET_TYPEOF(s, LANGSXP);
	SETCAR(t, install("res")); t = CDR(t);
	SETCAR(t, raster);
	ans=eval(s, rho);
	UNPROTECT(1);
	return ans;
}

SEXP getRasterValues(SEXP raster,SEXP rho) {
	SEXP ans,s,t;
	if(!isEnvironment(rho)) error("'rho' should be an environment");

	t = s = PROTECT(allocList(2));
	SET_TYPEOF(s, LANGSXP);
	SETCAR(t, install("values")); t = CDR(t);
	SETCAR(t, raster);
	ans=eval(s, rho);
	UNPROTECT(1);
	return ans;
}

/*
* Import a raster and prepare it for fast use
* WARNING: this function intentionally leaves one SEXP protected
*/
RASTER *openRaster(SEXP raster, SEXP rho) {
	RASTER *out = malloc(sizeof(RASTER));
	SEXP dim, extent;
// this will not be unprotected, to avoid having to copy the values to a new
// variable
	out->pvalues = getRasterValues(raster, rho);
	R_PreserveObject(out->pvalues);
	out->values = NUMERIC_POINTER(out->pvalues);

	PROTECT(dim = getRasterDim(raster, rho));
	PROTECT(extent = getRasterExtent(raster, rho));
	out->nrows = NUMERIC_POINTER(dim)[0];
	out->ncols = NUMERIC_POINTER(dim)[1];
	out->xmin = NUMERIC_POINTER(GET_SLOT(extent, SCALARCHAR("xmin")))[0];
	out->ymin = NUMERIC_POINTER(GET_SLOT(extent, SCALARCHAR("ymin")))[0];
	out->xmax = NUMERIC_POINTER(GET_SLOT(extent, SCALARCHAR("xmax")))[0];
	out->ymax = NUMERIC_POINTER(GET_SLOT(extent, SCALARCHAR("ymax")))[0];
	out->width = out->xmax - out->xmin;
	out->height = out->ymax - out->ymin;
	out->xscale = (float) out->ncols / out->width;
	out->yscale = (float) out->nrows / out->height;
	out->ncells = LENGTH(out->pvalues);
	
	UNPROTECT(2);
	return out;
}

void closeRaster(RASTER *raster) {
	R_ReleaseObject(raster->pvalues);
	free(raster);
}

/*
* Fast extraction of raster values at a given coordinate
*/
inline double extractRasterValue(const RASTER *raster, float x, float y) {
	if(x < raster->xmin || y < raster->ymin || x >= raster->xmax
		|| y >= raster->ymax) return NA_REAL;
	int index = (int) ((raster->ymax - y) * raster->yscale) * raster->ncols
		+ (int) ((x - raster->xmin) * raster->xscale);
	if(index >= raster->ncells) return NA_REAL;
	return raster->values[index];
}

/*
* Fast extraction of raster values at a given coordinate
* guaranteed not to return NaN
*/
inline double extractRasterValueNoNaN(const RASTER *raster,float x,float y) {
	if(x < raster->xmin || y < raster->ymin || x >= raster->xmax
		|| y >= raster->ymax) return 1;
	int index = (int) ((raster->ymax - y) * raster->yscale) * raster->ncols
		+ (int) ((x - raster->xmin) * raster->xscale);
	if(index >= raster->ncells) return 1;
	double tmp = raster->values[index];
	return isnan(tmp) ? 1 : tmp;
}

