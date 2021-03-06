#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h> // for NULL  // Probably unnecessary

#include "mpl.h"
#include "RMorphyUtils.h"
#include "RMorphy.h"

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
  R_RegisterCCallable("inapplicable", "MORPHYLENGTH", (DL_FUNC) &MORPHYLENGTH);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
