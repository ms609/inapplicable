#include <stdlib.h>
#include "RMorphy.h"

// Random number generator from http://www.cse.yorku.ca/~oz/marsaglia-rng.html
// 1+MWC%10 generates an integer from 1 to 10
#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define random_int ((znew<<16)+wnew)

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

// parent_of, left and right have been initialized with a two-taxon tree with tips 0 & 1
// left and right point n_tip _before_ left and right, so we don't need to subtract n_tip each time
// We arbitrarily choose to root our tree on tip 0, so never add to that edge or the 
// "dummy" root edge.
void build_tree(int *parent_of, int *left, int *right, const int *n_tip) {
  int i, addition_point, new_node;
  for (i = 3; i < *n_tip; i++) {
    new_node = i + *n_tip - 1;
    addition_point = 1 + random_int % (i + i - 3); // +1 to avoid edge 0
    if (addition_point < i) { // Adding below a tip
      insert_in_order(parent_of, left, right, &addition_point, &new_node, &i);
    } else { // Adding below an existing node
      addition_point += *n_tip - i + 1; // + 1 to avoid dummy root
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
  int  *replacement_array = malloc((*n_tip - 1)           * sizeof(int)),
              *parent_ref = malloc((*n_tip + *n_tip - 1)  * sizeof(int)),
              *left_array = malloc((*n_tip - 1)           * sizeof(int)),
             *right_array = malloc((*n_tip - 1)           * sizeof(int)),
      *replacement_number = replacement_array - *n_tip,
                *left_ref = left_array        - *n_tip,
               *right_ref = right_array       - *n_tip,
                        i = *n_tip,
                        j = *n_tip + 1;
  replacement_number[*n_tip] = *n_tip;
  move_to_node(&i, parent_of, left, right, replacement_number, &j, n_tip);
  
  for (i = 0; i < *n_tip; i++) {
    parent_ref[i] = parent_of[i];
  }
  for (i = *n_tip; i < (*n_tip + *n_tip - 1); i++) {
    parent_ref[i] = parent_of[i];
    left_ref  [i] = left[i];
    right_ref [i] = right[i];
  }
  
  for (i = 0; i < *n_tip; i++) {
    // Tips have not been renumbered; they are easy
    parent_of[i] = replacement_number[parent_ref[i]];
  }
  for (i = *n_tip; i < (*n_tip + *n_tip - 1); i++) {
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
  Morphy handl = R_ExternalPtrAddr(MorphyHandl);
  SEXP RESULT = PROTECT(allocVector(INTSXP, 1));
  int *score;
  score = INTEGER(RESULT);
  *score = 0;
  if (n_tip < 2) {
    INTEGER(RESULT)[0] = 0;
    UNPROTECT(1);
    return(RESULT);
  }
  
  int *parent_of = malloc(n_tip + n_tip - 1 * sizeof(int)),
           *left = malloc(n_tip - 1         * sizeof(int)),
          *right = malloc(n_tip - 1         * sizeof(int));
          
  if (n_tip < 3) {
        // Initialize with 2-tip tree
      parent_of[0] = n_tip;
      parent_of[1] = n_tip;
  parent_of[n_tip] = n_tip; // Root is its own parent
           left[0] = 0;
          right[0] = 1;
  } else {
    // Initialize with 3-tip tree, arbitrarily rooted on tip 0
        parent_of[0] = n_tip;
        parent_of[1] = n_tip + 1;
        parent_of[2] = n_tip + 1;
    parent_of[n_tip] = n_tip; // Root is its own parent
             left[0] = 0;
             left[1] = 1;
            right[0] = n_tip + 1;
            right[1] = 2;
  }
  if (n_tip > 3) {    
    build_tree(parent_of, left - n_tip, right - n_tip, &n_tip);
    renumber_postorder(parent_of, left - n_tip, right - n_tip, &n_tip);
  }
  
  morphy_length(parent_of, left, right, handl, score); 
  
  free(parent_of);
  free(right);
  free(left);
  UNPROTECT(1);
  return(RESULT);
}
