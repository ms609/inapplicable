#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "mpl.h"
#include "RMorphyUtils.h"



SEXP MORPHYLENGTH(SEXP R_ancestors, SEXP R_left, SEXP R_right, SEXP MorphyHandl) {
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  const int n_taxa = mpl_get_numtaxa(handl); 
  const int n_internal = mpl_get_num_internal_nodes(handl);
  const int root_node = n_taxa;
  const int max_node = n_taxa + n_internal - 1;
  
  // R_descendants and R_ancestors have already had one subtracted to convert them to an index 
  const int *ancestor=INTEGER(R_ancestors), *left=INTEGER(R_descendants), 
            *right=INTEGER(R_right); // INTEGER gives pointer to first element of length n R vector
  
  // Declare and protect result, to return to R
  SEXP Rres = PROTECT(allocVector(INTSXP, 1));
  
  // Initialize return variables
  int i;
  int *pscore_temp;
  pscore_temp = INTEGER(Rres);
  *pscore_temp = 0;
   
  for (i = max_node - 1; i > n_taxa; i--) { // First Downpass 
    *pscore_temp += mpl_first_down_recon(i,  left[i - n_taxa], right[i - n_taxa], handl);
  }
  mpl_update_lower_root(max_node, root_node, handl); // The index of the lower dummy root node = max_node/
                                                     // TODO we can alternatively pass the root node as its own ancestor.
 
  for (i = root_node; i < max_node; i++) { // First uppass: internal nodes
    *pscore_temp += mpl_first_up_recon(i, left[i - n_taxa], right[i - n_taxa], ancestor[i], handl);
  }
  for (i = 0; i < n_taxa; i++) { // First uppass: update tips
    mpl_update_tip(i, ancestor[i], handl);
  }
  
  for (i = max_node - 1; i > n_taxa; i--) { // Second Downpass 
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