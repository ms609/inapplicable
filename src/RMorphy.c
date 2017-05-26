// Useful resource for R-C interface: http://adv-r.had.co.nz/C-interface.html
// Useful tool for debugging: function "cfunction" in R package "inline"

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "mpl.h"
#include "RMorphyUtils.h"


// Morphy handler for R
// Simplest implementation (good place to start!): 
// Return the number of steps associated with each character 
// This can then be multiplied by a vector of character weights in R
// This will probably be quick enough to allow tree search to proceed in R, even without 
// subsequent optimisations (e.g. using 'normal' fitch optimization for characters with <2 inapp tokens)

SEXP RMorphy(SEXP r_n_char, SEXP r_n_taxa, SEXP r_descendants, SEXP r_ancestors, SEXP r_rawmatrix) {
  
  // Convert arguments passed from R into C variables
  
  int n_char=asInteger(r_n_char), n_taxa=asInteger(r_n_taxa);  // asInteger converts length 1 R vector into scalar
  // r_descendants and r_ancestors have already had one subtracted to convert them to an index 
  int *descendants=INTEGER(r_descendants), *ancestors=INTEGER(r_ancestors);  // INTEGER gives pointer to first element of length n R vector
  const char *rawmatrix=CHAR(asChar(r_rawmatrix));
 
  Rprintf(" - character string = %s\n", rawmatrix);
  
  // Calculate relevant properties of the tree and dataset
  int n_internal = n_taxa; // One more than you might expect because there's a dummy root node
  int root_node = n_taxa;
  int max_node = n_taxa + n_internal - 1L;
  Rprintf(" - %i nodes + dummy root, of which %i are internal (plus dummy root), %i are tips\n", max_node, n_internal - 1L, n_taxa);
  
  // Declare and protect result, to return to R
  SEXP RESULT, pscore, node_children;
  PROTECT(RESULT = allocVector(VECSXP, 2));
  PROTECT(pscore = allocVector(INTSXP, 1));
  PROTECT(node_children = allocVector(INTSXP, 2));
  
  // Initialize return variables
  int i;
  int *pscore_temp, *children_temp;
  pscore_temp = INTEGER(pscore);
  *pscore_temp = 0;
  children_temp = INTEGER(node_children);
  for (i = 0; i < 2; i++) {
    children_temp[i] = 0;
  }
    
  // This is just to test that the ancestor identification is working as expected
  int root_node_left_child = n_taxa + 1;
  children_temp[0] = descendants[(root_node_left_child - n_taxa) * 2];
  children_temp[1] = descendants[((root_node_left_child - n_taxa) * 2)+ 1];
  
  Rprintf(" - Root node's left child has nodes %i and %i as its left and right children\n",
    children_temp[0], children_temp[1]
  );
  
  Morphy handl = mpl_new_Morphy();
  mpl_init_Morphy(n_taxa, n_char, handl);
 
  for (i = 0; i < n_char; i++) {
    mpl_set_parsim_t(i, FITCH_T, handl);
  }
  mpl_set_num_internal_nodes(n_internal, handl);
  mpl_attach_rawdata(rawmatrix, handl);
  mpl_apply_tipdata(handl);
  Rprintf(" - Initialized morphy object with %i characters and %i taxa.\n - Begin traversal:\n", mpl_get_num_charac(handl), mpl_get_numtaxa(handl));
  
 
  for (i = max_node - 1L; i > n_taxa; i--) { // First Downpass 
  // First node is called 11 in R, it's numbered 10 in C
  // Max node is 11, so we knock one off to give the node's index, i = 10
  // Stop before i == n_taxa (which would be the root node)
  
  // i == 11
  // If i == 6 (root node), descendants will be in positions [0, 1].  i - n_taxa = 0
  // If i == 11, descendants will be in positions [8, 9]. i - n_taxa = 4
  
    Rprintf("   - Reconstructing  node %i < %i,%i.\n", i, descendants[(i - n_taxa) * 2], descendants[((i - n_taxa) * 2) + 1]);
    *pscore_temp += mpl_first_down_recon(i,  descendants[(i - n_taxa) * 2], descendants[((i - n_taxa) * 2) + 1], handl);
  }
  mpl_update_lower_root(max_node, root_node, handl);
 
  for (i = root_node; i < max_node; i++) { // First uppass: internal nodes
    *pscore_temp += mpl_first_up_recon(i, descendants[(i - n_taxa) * 2], descendants[((i - n_taxa) * 2) + 1], ancestors[i], handl);
  }
  for (i = 0; i < n_taxa; i++) { // First uppass: update tips
    mpl_update_tip(i, ancestors[i], handl);
  }
  
  for (i = max_node - 1L; i > n_taxa; i--) { // Second Downpass 
    *pscore_temp += mpl_second_down_recon(i, descendants[(i - n_taxa) * 2], descendants[((i - n_taxa) * 2) + 1], handl);
  }
 
  for (i = n_taxa; i < max_node; i++) { // Second uppass: internal nodes
    *pscore_temp += mpl_second_up_recon(i, descendants[(i - n_taxa) * 2], descendants[((i - n_taxa) * 2) + 1], ancestors[i], handl);
  }
  for (i = 0; i < n_taxa; i++) { // Second uppass: finalize tips (fwiw)
    mpl_finalize_tip(i, ancestors[i], handl);
  }
  
  Rprintf("\n Traversal complete.\n - Deleting Morphy handle\n");
  mpl_delete_Morphy(handl);
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, node_children);
  UNPROTECT(3);
  return (RESULT);
  
  return R_NilValue;
}