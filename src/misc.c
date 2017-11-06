#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

/**
* Get the list element named str, or return NULL
*/
SEXP getListElement(SEXP list, const char *str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	for (int i = 0; i < length(list); i++)
	 if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
		elmt = VECTOR_ELT(list, i);
		break;
	 }
	return elmt;
}

