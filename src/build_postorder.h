#include <stdlib.h>
void insert_in_order (int *parent_of, int *left, int *right, 
                      const int *addition_point, const int *new_node
                      const int *new_tip) {
  const int old_parent = parent_of[*addition_point];
  if (left[old_parent] == *addition_point) {
    left[old_parent] = *new_node;
  } else {
    // The same, but on the right
    right[old_parent] = *new_node;
  }
  left[*new_node] = *new_tip;
  parent_of[*new_tip] = *new_node;
  right[*new_node] = *addition_point;
  parent_of[*addition_point] = *new_node;
}

void swap_nodes(const int node1, const int *node2, parent_of, left, right) {
  // TODO WRITE THIS FUNCTION
}

static R_NativePrimitiveArgType build_postorder_tree_t[] = {
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
// parent_of, left and right have been initialized with a two-taxon tree with tips 0 & 1
// left and right point n_tip _before_ left and right, so we don't need to subtract n_tip each time
extern void build_postorder_tree(int *parent_of, int *left, int *right, const int *n_tip)
{
  int i, addition_point, new_node;
  for (i = 2; i < *n_tip; i++) {
    new_node = i + *n_tip;
    addition_point = rand() % (i + i - 1L);
    if (addition_point < *n_tip) { // Adding below a tip
      insert_in_order(parent_of, left, right, &addition_point, &new_node, &i);
    } else { // Adding below an existing node
      addition_point += n_tip;
      insert_in_order(parent_of, left, right, &addition_point, &new_node, &i);
      if (parent_of[addition_point] < new_node) {
        swap_nodes(parent_of[addition_point], new_node, parent_of, left, right);
      }
    }
  }
}

extern SEXP BUILD_POSTORDER(SEXP ntip, SEXP MorphyHandl) {
  // tipnames run from 0 to nTip - 1, in random order
  const int n_tip = INTEGER(ntip)[0];
  int i, score;
  PROTECT(RESULT = allocVector(INTSXP, 1L));
  if (n_tip < 2) {
    INTEGER(RESULT)[0] = 0;
    UNPROTECT(1);
    return(RESULT);
  }
  
  int *parent_of = calloc(n_tip + n_tip - 1L, sizeof(int)),
           *left = calloc(n_tip - 1L, sizeof(int)),
          *right = calloc(n_tip - 1L, sizeof(int));
  // Initialize with 
      parent_of[0] = n_tip;
      parent_of[1] = n_tip;
  parent_of[n_tip] = n_tip; // Root is its own parent
          right[0] = 0; // Already true
           left[0] = 1;
  
  build_postorder_tree(parent_of, left - n_tip, right - n_tip, &n_tip);
  // Then score:
  // score = MORPHYLENGTH(parent_of, left, right, MorphyHandl); 
  // Can you send R objects as SEXPs?
  // If not, we'll need to hollow out MORPHYLENGTH to create a callable.
  
  INTEGER(RESULT)[0] score;
  UNPROTECT(1);
  return(RESULT);
}
