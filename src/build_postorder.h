#include <stdlib.h>


void insert_in_order (int *parent_of, int *left, int *right, 
                      const int *addition_point, const int *new_node,
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
  parent_of[*new_node] = old_parent;
}

void insert_and_reorder (int *parent_of, int *left, int *right, 
                      const int *addition_point, const int *new_node,
                      const int *new_tip) {
   left[*new_node] = left [*addition_point];
  right[*new_node] = right[*addition_point];
  parent_of[left [*new_node]] = *new_node;
  parent_of[right[*new_node]] = *new_node;
  left[*addition_point] = *new_tip;
  parent_of[*new_tip] = *addition_point;
  parent_of[*new_node] = *addition_point;
  right[*addition_point] = *new_node;
}

// parent_of, left and right have been initialized with a two-taxon tree with tips 0 & 1
// left and right point n_tip _before_ left and right, so we don't need to subtract n_tip each time
void build_tree(int *parent_of, int *left, int *right, const int *n_tip) {
  int i, addition_point, new_node;
  for (i = 2; i < *n_tip; i++) {
    new_node = i + *n_tip - 1L;
    addition_point = rand() % (i + i - 1L);
    if (addition_point < i) { // Adding below a tip
      insert_in_order(parent_of, left, right, &addition_point, &new_node, &i);
    } else if (addition_point == i) { // Adding below root node
      insert_and_reorder(parent_of, left, right, n_tip, &new_node, &i);
    } else { // Adding below an existing node
      addition_point += *n_tip - i;
      insert_in_order(parent_of, left, right, &addition_point, &new_node, &i);
    }
  }
}

void move_to_node(const int *node, const int *parent_of, const int *left, const int *right, 
                  int *replacement, int *next_label, const int *n_tip) {
  if (left [*node] > *n_tip) { // won't be equal, as that the root is no-one's descendant
    replacement[left[*node]] = (*next_label)++;
    move_to_node(&(left[*node]), parent_of, left, right, replacement, next_label, n_tip);
  }
  if (right[*node] > *n_tip) {
    replacement[right[*node]] = (*next_label)++;
    move_to_node(&right[*node], parent_of, left, right, replacement, next_label, n_tip);    
  }
}

void renumber_postorder(int *parent_of, int *left, int *right, const int *n_tip) {
  int  *replacement_array = malloc((*n_tip - 1L)           * sizeof(int)),
              *parent_ref = malloc((*n_tip + *n_tip - 1L)  * sizeof(int)),
              *left_array = malloc((*n_tip - 1L)           * sizeof(int)),
             *right_array = malloc((*n_tip - 1L)           * sizeof(int)),
      *replacement_number = replacement_array - *n_tip,
                *left_ref = left_array        - *n_tip,
               *right_ref = right_array       - *n_tip,
                        i = *n_tip,
                        j = *n_tip + 1L;
  replacement_number[*n_tip] = *n_tip;
  move_to_node(&i, parent_of, left, right, replacement_number, &j, n_tip);
  
  for (i = 0; i < *n_tip; i++) {
    parent_ref[i] = parent_of[i];
  }
  for (i = *n_tip; i < (*n_tip + *n_tip - 1L); i++) {
    parent_ref[i] = parent_of[i];
    left_ref  [i] = left[i];
    right_ref [i] = right[i];
  }
  
  for (i = 0; i < *n_tip; i++) {
    // Tips have not been renumbered; they are easy
    parent_of[i] = replacement_number[parent_ref[i]];
  }
  for (i = *n_tip; i < (*n_tip + *n_tip - 1L); i++) {
    // Nodes may have been renumbered; make sure we use the new numbers.
    parent_of[i] = replacement_number[parent_ref[replacement_number[i]]];
    left[i] = (left_ref[replacement_number[i]] > *n_tip) ?
                replacement_number[left_ref[replacement_number[i]]] :
                left_ref[replacement_number[i]];
    right[i] = (right_ref[replacement_number[i]] > *n_tip) ?
                replacement_number[right_ref[replacement_number[i]]] :
                right_ref[replacement_number[i]];
  }
  free(replacement_array);
  free(right_array);
  free(left_array);
  free(parent_ref);
}

extern SEXP BUILD_POSTORDER(SEXP ntip, SEXP MorphyHandl) {
  // tipnames run from 0 to nTip - 1, in random order
  const int n_tip = INTEGER(ntip)[0];
  //int i, score;
  //PROTECT(RESULT = allocVector(INTSXP, 1L));
  SEXP RESULT, PARENT_OF, RIGHT, LEFT;
  PROTECT(RESULT = allocVector(VECSXP, 3L));
  PROTECT(PARENT_OF = allocVector(INTSXP, n_tip + n_tip - 1L));
  PROTECT(LEFT      = allocVector(INTSXP, n_tip - 1L));
  PROTECT(RIGHT     = allocVector(INTSXP, n_tip - 1L));
  if (n_tip < 2) {
    INTEGER(RESULT)[0] = 0;
    UNPROTECT(1);
    return(RESULT);
  }
  
  // int *parent_of = calloc(n_tip + n_tip - 1L, sizeof(int)),
  //          *left = calloc(n_tip - 1L, sizeof(int)),
  //         *right = calloc(n_tip - 1L, sizeof(int));
  
  int *parent_of = INTEGER(PARENT_OF),
          *right = INTEGER(RIGHT),
           *left = INTEGER(LEFT);
  // Initialize with 
      parent_of[0] = n_tip;
      parent_of[1] = n_tip;
  parent_of[n_tip] = n_tip; // Root is its own parent
          right[0] = 0; // Already true
           left[0] = 1;
  
  build_tree(parent_of, left - n_tip, right - n_tip, &n_tip);
  renumber_postorder(parent_of, left - n_tip, right - n_tip, &n_tip);
  
  // score = MORPHYLENGTH(parent_of, left, right, MorphyHandl); 
  // Can you send R objects as SEXPs?
  // If not, we'll need to hollow out MORPHYLENGTH to create a callable.
  
  // free (*parent_of);
  // free (*right);
  // free (*left);
  
  SET_VECTOR_ELT(RESULT, 0, PARENT_OF);
  SET_VECTOR_ELT(RESULT, 1, LEFT);
  SET_VECTOR_ELT(RESULT, 2, RIGHT);
  //INTEGER(RESULT)[0] = score;
  //UNPROTECT(1);
  UNPROTECT(4);
  return(RESULT);
}
