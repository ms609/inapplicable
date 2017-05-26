#include <R.h>
#include <Rinternals.h>
#include "mpl.h"
#include "RMorphyUtils.h"

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
