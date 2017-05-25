#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


// Morphy handler for R
SEXP morphy_handler_R(SEXP matrix_R, SEXP ntaxa_R, SEXP nchar_R, SEXP list_anc, SEXP list_child, SEXP weight) {

    // Input variables
    int ntaxa = asInteger(ntaxa_R);
    int nchar = asInteger(nchar_R);
    const char *matrix = CHAR(asChar(matrix_R)); // As the matrix will be sent in the format of a single string

    // Vectors from R
    
    list_anc = coerceVector(list_anc, INTSXP);
    PROTECT(list_anc);
    list_child = coerceVector(list_child, INTSXP);
    PROTECT(list_child);
    weight = coerceVector(weight, INTSXP);
    PROTECT(weight);

    // Variables from C
    int nedges = 2 * ntaxa - 2;
    int nnodes = ntaxa - 1;
    char* symbols = NULL;
    int i = 0;

    // Create a morphy handle Here


    

    // Do some stuff on the handle
    Rprintf("ntaxa  = %i \n", ntaxa);
    Rprintf("nchar  = %i \n", nchar);

    printf("matrix  = ");
    for(i = 0; i < ntaxa*nchar; ++i) {
        Rprintf("%c", CHAR(asChar(matrix))[i]) ;
    }
    Rprintf("\n");

    Rprintf("ancestors  = ");
    for(i = 0; i < nnodes+ntaxa; ++i) {
        Rprintf("%i", INTEGER(list_anc)[i]) ;
    }
    Rprintf("\n");

    Rprintf("children  = ");
    for(i = 0; i < nnodes*2; ++i) {
        Rprintf("%i", INTEGER(list_child)[i]) ;
    }
    Rprintf("\n");

    Rprintf("weights  = ");
    for(i = 0; i < nchar; ++i) {
        Rprintf("%i", INTEGER(weight)[i]) ;
    }
    Rprintf("\n");

    Rprintf("nedges  = %i \n", nedges);
    Rprintf("nnodes  = %i \n", nnodes);


    // Reconvert outputs


    // Some silly output
    SEXP output = PROTECT(allocVector(REALSXP, 1));

    UNPROTECT(4);

    // Return
    return output;
}