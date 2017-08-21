#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "mpl.h"
#include "RMorphyUtils.h"
#include <R_ext/Rdynload.h>

SEXP _R_wrap_mpl_new_Morphy(void)
{
    Morphy new = mpl_new_Morphy();
    SEXP result = R_MakeExternalPtr(new, R_NilValue, R_NilValue);
    
    return result;
}

SEXP _R_wrap_mpl_delete_Morphy(SEXP MorphyHandl)
{
    int ret = 0;
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    Morphy handl = R_ExternalPtrAddr(MorphyHandl);
    ret = mpl_delete_Morphy(handl);
    INTEGER(Rret)[0] = ret;
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_init_Morphy(SEXP Rntax, SEXP Rnchar, SEXP MorphHandl)
{
    int ret = 0;
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    Morphy handl = R_ExternalPtrAddr(MorphHandl);
    int ntax = INTEGER(Rntax)[0];
    int nchar = INTEGER(Rnchar)[0];

    ret = mpl_init_Morphy(ntax, nchar, handl);

    INTEGER(Rret)[0] = ret;
    UNPROTECT(1);
    return Rret;
}
    
SEXP _R_wrap_mpl_get_numtaxa(SEXP MorphHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_get_numtaxa(R_ExternalPtrAddr(MorphHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_num_charac(SEXP MorphHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    
    INTEGER(Rret)[0] = mpl_get_num_charac(R_ExternalPtrAddr(MorphHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_attach_symbols(SEXP Rsymbols, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    int Mret = 0;
    const char *Msymbols = CHAR(asChar(Rsymbols));

    Mret = mpl_attach_symbols(Msymbols, R_ExternalPtrAddr(MorphyHandl));
    INTEGER(Rret)[0] = Mret;

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_symbols(SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(STRSXP, 1));

    char* symbols = mpl_get_symbols(R_ExternalPtrAddr(MorphyHandl));

    Rret = mkString(symbols);
    
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_attach_rawdata(SEXP Rmatrix, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    int Mret = 0;
    const char *Mmatrix = CHAR(asChar(Rmatrix));

    Mret = mpl_attach_rawdata(Mmatrix, R_ExternalPtrAddr(MorphyHandl));
    
    INTEGER(Rret)[0] = Mret;
    UNPROTECT(1);

    return Rret;
}

SEXP _R_wrap_mpl_delete_rawdata(SEXP MorphyHandl)
{
	SEXP Rret = PROTECT(allocVector(INTSXP, 1));
	Morphy handl = R_ExternalPtrAddr(MorphyHandl);
    int ret = 0;

    ret = mpl_delete_rawdata(handl);
    
    INTEGER(Rret)[0] = ret;
    UNPROTECT(1);

    return Rret;
}

SEXP _R_wrap_mpl_set_parsim_t(SEXP RcharID, SEXP Rchtype, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));
    MPLchtype chtype;
    int Mret = 0;

    const char* chtypename = CHAR(asChar(Rchtype));
    
    chtype = _R_mpl_str2chtype(chtypename);
    Mret = mpl_set_parsim_t(INTEGER(RcharID)[0], chtype,
                            R_ExternalPtrAddr(MorphyHandl));

    INTEGER(Rret)[0] = Mret;

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_set_gaphandl(SEXP Rgaptype, SEXP MorphyHandl)
{
	SEXP Rret = PROTECT(allocVector(INTSXP, 1));
	MPLgap_t gaptype;
	int Mret = 0;

	const char* gaptypename = CHAR(asChar(Rgaptype));
	
	gaptype = _R_mpl_str2gaptype(gaptypename);
	Mret = mpl_set_gaphandl(gaptype, R_ExternalPtrAddr(MorphyHandl));

	INTEGER(Rret)[0] = Mret;

	UNPROTECT(1);
	return Rret;
}

SEXP _R_wrap_mpl_set_num_internal_nodes(SEXP Rnnodes, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_set_num_internal_nodes
                        (INTEGER(Rnnodes)[0], 
                         R_ExternalPtrAddr(MorphyHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_get_num_internal_nodes(SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_get_num_internal_nodes
                        (R_ExternalPtrAddr(MorphyHandl));

    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_apply_tipdata(SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = mpl_apply_tipdata(R_ExternalPtrAddr(MorphyHandl));
    
    UNPROTECT(1);    
    return Rret;
}


SEXP _R_wrap_mpl_set_charac_weight(SEXP RcharID, SEXP Rweight, SEXP MorphyHandl)
{
  SEXP Rret = PROTECT(allocVector(INTSXP, 1));
  INTEGER(Rret)[0] = mpl_set_charac_weight(INTEGER(RcharID)[0], REAL(Rweight)[0],
                                           R_ExternalPtrAddr(MorphyHandl));
  UNPROTECT(1);
  return Rret;
}

SEXP _R_wrap_mpl_get_charac_weight(SEXP RcharID, SEXP MorphyHandl) 
{
  SEXP Rret, Wapprox, Wexact;
  PROTECT(Rret = allocVector(VECSXP, 2));
  PROTECT(Wapprox = allocVector(INTSXP, 1));
  PROTECT(Wexact  = allocVector(REALSXP, 1));
  
  INTEGER(Wapprox)[0] = mpl_get_charac_weight(REAL(Wexact), INTEGER(RcharID)[0],
                                              R_ExternalPtrAddr(MorphyHandl));
  SET_VECTOR_ELT(Rret, 0, Wapprox);
  SET_VECTOR_ELT(Rret, 1, Wexact);
  UNPROTECT(3);
  return Rret;
}

SEXP _R_wrap_mpl_first_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_first_down_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                         INTEGER(Rright_id)[0],
                         R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_first_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP Ranc_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_first_up_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       INTEGER(Rright_id)[0], INTEGER(Ranc_id)[0],
                       R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_second_down_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_second_down_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       	  INTEGER(Rright_id)[0],
                          R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_second_up_recon(SEXP Rnode_id, SEXP Rleft_id, SEXP Rright_id, SEXP Ranc_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_second_up_recon(INTEGER(Rnode_id)[0], INTEGER(Rleft_id)[0],
                       	INTEGER(Rright_id)[0], INTEGER(Ranc_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_update_tip(SEXP tip_id, SEXP anc_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_update_tip(INTEGER(tip_id)[0], INTEGER(anc_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP _R_wrap_mpl_update_lower_root(SEXP lower_id, SEXP upper_id, SEXP MorphyHandl)
{
    SEXP Rret = PROTECT(allocVector(INTSXP, 1));

    INTEGER(Rret)[0] = 
    mpl_update_lower_root(INTEGER(lower_id)[0], INTEGER(upper_id)[0],
                        R_ExternalPtrAddr(MorphyHandl));
    UNPROTECT(1);
    return Rret;
}

SEXP MORPHYLENGTH(SEXP R_ancestors, SEXP R_left, SEXP R_right, SEXP MorphyHandl) {
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  const int n_taxa = mpl_get_numtaxa(handl); 
  const int n_internal = mpl_get_num_internal_nodes(handl);
  const int root_node = n_taxa;
  const int max_node = n_taxa + n_internal;
  
  // R_descendants and R_ancestors have already had one subtracted to convert them to an index 
  const int *ancestor=INTEGER(R_ancestors), *left=INTEGER(R_left), 
            *right=INTEGER(R_right); // INTEGER gives pointer to first element of length n R vector
  
  // Declare and protect result, to return to R
  SEXP Rres = PROTECT(allocVector(INTSXP, 1));
  
  // Initialize return variables
  int i;
  int *pscore_temp;
  pscore_temp = INTEGER(Rres);
  *pscore_temp = 0;
   
  for (i = max_node - 1; i >= n_taxa; i--) { // First Downpass 
    *pscore_temp += mpl_first_down_recon(i,  left[i - n_taxa], right[i - n_taxa], handl);
    //Rprintf("Downpass on node %i -< %i,%i ... pscore is %i\n", i, left[i-n_taxa], right[i-n_taxa], *pscore_temp);
    
  }
  mpl_update_lower_root(root_node, root_node, handl); // We could use a spare internal node with index = max_node as a dummy root node.
                                                      // Or we can just pass the root node as its own ancestor.
 
  for (i = root_node; i < max_node; i++) { // First uppass: internal nodes
    *pscore_temp += mpl_first_up_recon(i, left[i - n_taxa], right[i - n_taxa], ancestor[i], handl);
    //Rprintf("Uppass on node %i -< %i,%i ... pscore is %i\n", i, left[i-n_taxa], right[i-n_taxa], *pscore_temp);
  }
  for (i = 0; i < n_taxa; i++) { // First uppass: update tips
    mpl_update_tip(i, ancestor[i], handl);
  }
  
  for (i = max_node - 1; i >= n_taxa; i--) { // Second Downpass 
    *pscore_temp += mpl_second_down_recon(i, left[i - n_taxa], right[i - n_taxa], handl);
  }
 
  for (i = n_taxa; i < max_node; i++) { // Second uppass: internal nodes
    *pscore_temp += mpl_second_up_recon(i, left[i - n_taxa], right[i - n_taxa], ancestor[i], handl);
  }
  for (i = 0; i < n_taxa; i++) { // Second uppass: finalize tips (fwiw)
    mpl_finalize_tip(i, ancestor[i], handl);
  }
  UNPROTECT(1);
  return Rres;
}


static const R_CallMethodDef callMethods[] = {
  {"_R_wrap_mpl_new_Morphy",        (DL_FUNC) &_R_wrap_mpl_new_Morphy, 0},
  {"_R_wrap_mpl_delete_Morphy",     (DL_FUNC) &_R_wrap_mpl_delete_Morphy, 1},
  {"_R_wrap_mpl_init_Morphy",       (DL_FUNC) &_R_wrap_mpl_init_Morphy, 3},
  {"_R_wrap_mpl_get_numtaxa",       (DL_FUNC) &_R_wrap_mpl_get_numtaxa, 1},
  {"_R_wrap_mpl_get_num_charac",    (DL_FUNC) &_R_wrap_mpl_get_num_charac, 1},
  {"_R_wrap_mpl_attach_symbols",    (DL_FUNC) &_R_wrap_mpl_attach_symbols, 2},
  {"_R_wrap_mpl_get_symbols",       (DL_FUNC) &_R_wrap_mpl_get_symbols, 1},
  {"_R_wrap_mpl_attach_rawdata",    (DL_FUNC) &_R_wrap_mpl_attach_rawdata, 2},
  {"_R_wrap_mpl_delete_rawdata",    (DL_FUNC) &_R_wrap_mpl_delete_rawdata, 1},
  {"_R_wrap_mpl_set_parsim_t",      (DL_FUNC) &_R_wrap_mpl_set_parsim_t, 3},
  {"_R_wrap_mpl_set_gaphandl",      (DL_FUNC) &_R_wrap_mpl_set_gaphandl, 2},
  {"_R_wrap_mpl_set_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_set_num_internal_nodes, 2},
  {"_R_wrap_mpl_get_num_internal_nodes", (DL_FUNC) &_R_wrap_mpl_get_num_internal_nodes, 1},
  {"_R_wrap_mpl_apply_tipdata",     (DL_FUNC) &_R_wrap_mpl_apply_tipdata, 1},
  {"_R_wrap_mpl_set_charac_weight", (DL_FUNC) &_R_wrap_mpl_set_charac_weight, 3},
  {"_R_wrap_mpl_get_charac_weight", (DL_FUNC) &_R_wrap_mpl_get_charac_weight, 2},
  {"_R_wrap_mpl_first_down_recon",  (DL_FUNC) &_R_wrap_mpl_first_down_recon, 4},
  {"_R_wrap_mpl_first_up_recon",    (DL_FUNC) &_R_wrap_mpl_first_up_recon, 5},
  {"_R_wrap_mpl_second_down_recon", (DL_FUNC) &_R_wrap_mpl_second_down_recon, 4},
  {"_R_wrap_mpl_second_up_recon",   (DL_FUNC) &_R_wrap_mpl_second_up_recon, 5},
  {"_R_wrap_mpl_update_tip",        (DL_FUNC) &_R_wrap_mpl_update_tip, 3},
  {"_R_wrap_mpl_update_lower_root", (DL_FUNC) &_R_wrap_mpl_update_lower_root, 3},
  {"MORPHYLENGTH",                  (DL_FUNC) &MORPHYLENGTH, 4},
  {NULL, NULL, 0}
};

void R_init_inapplicable(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
 // R_RegisterCCallable("inapplicable", "MORPHYLENGTH", (DL_FUNC) &MORPHYLENGTH);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
