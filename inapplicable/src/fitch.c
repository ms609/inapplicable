#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void fitch_downnode_ia(int *dat1, int *dat2, int *n_rows, int *inapp) {
  int k, tmp;
  for (k = 0; k < (*n_rows); k++) {
    tmp = dat1[k] & dat2[k];
    if (tmp) {                                  // Children have tokens in Common
      if ((tmp == (*inapp))) {                  // Is - the only common token?
        if ((dat1[k] | dat2[k]) != (*inapp)) {  // Do either have an applicable token?
          if ((dat1[k] != (*inapp)) && (dat2[k] != (*inapp))) { // Do both have an applicable token?
            tmp = dat1[k] | dat2[k];            // Set tokens to those possible for left OR right child
          }
        }
      }
    } else {
      tmp = dat1[k] | dat2[k];                  // Set tokens to those possible for left OR right child
      if ((dat1[k] == *inapp) || (dat2[k] == *inapp)) { // Is one child's only possible token {-}?
      } else {
        if (tmp & (*inapp)) {                   // [Is {-} among this node's tokens?]
          tmp = tmp - (*inapp);                 // Delete {-} from this node's tokens
        }
      }
    }
    dat1[k] = tmp;
  }
}

void fitch_downpass_ia(int *dat, int *n_rows, int *parent, int *child, int *n_edge, int *inapp) {
  int i, ni, k;
  ni = 0;
  for (i = 0; i < *n_edge; i++) {
    if (ni == parent[i]){
      fitch_downnode_ia(&dat[(ni-1L) * (*n_rows)], &dat[(child[i]-1L) * (*n_rows)], n_rows, inapp);
    } else {
      ni = parent[i];
      for (k = 0; k < (*n_rows); k++) dat[(ni-1L)*(*n_rows) + k] = dat[(child[i]-1L)*(*n_rows) + k];                     
    }
  }
}

SEXP FITCHDOWNIA(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP n_edge, SEXP mx, SEXP q, SEXP inapp) {   
  int *data, *n_rows=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
  SEXP DAT;
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  data = INTEGER(DAT);
  for(i=0; i<(*n_rows * n); i++) data[i] = INTEGER(dat)[i];
  fitch_downpass_ia(data, n_rows, INTEGER(parent), INTEGER(child), INTEGER(n_edge), INTEGER(inapp));
  UNPROTECT(1);
  return(DAT);
}

void fitch_uproot_ia(int *this,                int *child_q, int *child_r, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, ancestor_k, q_inapp, r_inapp, common;
  for (k = 0; k < (*n_rows); k++) { // Next TS
    ancestor_k = this[k];
    if ((ancestor_k & (*inapp)) && (ancestor_k != (*inapp))) ancestor_k -= (*inapp);  // Remove {-} from parent node: Hennig's Auxiliary Principle
    r_inapp = (child_r[k] == *inapp);
    if ((ancestor_k != *inapp)                           // Parent is applicable
    && ((q_inapp = (child_q[k] == *inapp)) || r_inapp)) { // One child is inapplicable
      if (q_inapp && r_inapp) {
        child_q[k] = ancestor_k;
        child_r[k] = ancestor_k;
        this[k] = ancestor_k;
      } else if (q_inapp) { // Child Q is the inapplicable child
        if ((common = child_r[k] & ancestor_k)) { // Parent and applicable child have tokens in common
          child_q[k] = common;
          this[k] = common;
        } else {
          child_q[k] = child_r[k] | ancestor_k;
          if (child_q[k] & *inapp) child_q[k] -= *inapp;
          this[k] = child_q[k];          
          (pars[k])++;                              // Increase parsimony score by one
          (*w) += weight[k];                        // Increase parsimony score by one
        }
      } else { // child R is the inapplicable child
        if ((common = child_q[k] & ancestor_k)) { // Parent and applicable child have tokens in common
          child_r[k] = common;
          this[k] = common;
        } else {
          child_r[k] = child_q[k] | ancestor_k;
          if (child_r[k] & *inapp) child_r[k] -= *inapp;
          this[k] = child_r[k];
          (pars[k])++;                              // Increase parsimony score by one
          (*w) += weight[k];                        // Increase parsimony score by one
        }
      }
    } else if ((ancestor_k & this[k]) == ancestor_k) { // All parent node's tokens among this node's possible tokens
      this[k] = ancestor_k;                      // Set this node's tokens to parent's tokens
      if (((common = (child_q[k] & child_r[k])) && common != *inapp)   // Children have tokens in common, except {-}
        || child_q[k] == *inapp || child_r[k] == *inapp) // Either child == {-}
      {} else {
        if ((common = (child_q[k] | child_r[k])) & *inapp) this[k] = common - *inapp; else this[k] = common; // Tokens to union w/o -
        (pars[k])++;                              // Increase parsimony score by one
        (*w) += weight[k];                        // Increase parsimony score by one
      }
    } else if (child_q[k] & child_r[k]) {         // Children have tokens in common
      this[k] = this[k] | (ancestor_k & child_q[k]) | (ancestor_k & child_r[k]); 
                                                  // Add tokens common to parent and either child to this node
    } else {                                     
      this[k] = this[k] | ancestor_k;            // Add parent's tokens to this node's tokens
      if ((child_q[k] | child_r[k]) & *inapp) {} else {
        (pars[k])++;                              // Increase parsimony score by one
        (*w) += weight[k];                        // Increase parsimony score by one
      }
    }
  }
}

void fitch_upnode_ia(int *this, int *ancestor, int *child_q, int *child_r, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, q_inapp, r_inapp, common;
  for (k = 0; k < (*n_rows); k++) { // Next node in preorder
    r_inapp = (child_r[k] == *inapp);
    if ((ancestor[k] != *inapp)                           // Parent is applicable
    && ((q_inapp = (child_q[k] == *inapp)) || r_inapp)) { // One child is inapplicable
      if (q_inapp && r_inapp) {
        child_q[k] = ancestor[k];
        child_r[k] = ancestor[k];
        this[k] = ancestor[k];
      } else if (q_inapp) { // Child Q is the inapplicable child
        if ((common = child_r[k] & ancestor[k])) { // Parent and applicable child have tokens in common
          child_q[k] = common;
          this[k] = common;
        } else {
          child_q[k] = child_r[k] | ancestor[k];
          if (child_q[k] & *inapp) child_q[k] -= *inapp;
          this[k] = child_q[k];          
          (pars[k])++;                              // Increase parsimony score by one
          (*w) += weight[k];                        // Increase parsimony score by one
        }
      } else { // child R is the inapplicable child
        if ((common = child_q[k] & ancestor[k])) { // Parent and applicable child have tokens in common
          child_r[k] = common;
          this[k] = common;
        } else {
          child_r[k] = child_q[k] | ancestor[k];
          if (child_r[k] & *inapp) child_r[k] -= *inapp;
          this[k] = child_r[k];
          (pars[k])++;                              // Increase parsimony score by one
          (*w) += weight[k];                        // Increase parsimony score by one
        }
      }
    } else if ((ancestor[k] & this[k]) == ancestor[k]) { // All parent node's tokens among this node's possible tokens
      this[k] = ancestor[k];                      // Set this node's tokens to parent's tokens
      if (
      (
        !((common = (child_q[k] & child_r[k])) && common != *inapp)  // Children have tokens in common, except {-}: NO
      ) 
      && child_q[k] != *inapp && child_r[k] != *inapp
      ) { // Neither child == {-}: YES
        if ((common = (child_q[k] | child_r[k])) & *inapp) this[k] = common - *inapp; else this[k] = common; // Tokens to union w/o -
        (pars[k])++;                              // Increase parsimony score by one
        (*w) += weight[k];                        // Increase parsimony score by one
      }
    } else if (child_q[k] & child_r[k]) {         // Children have tokens in common
      this[k] = this[k] | (ancestor[k] & child_q[k]) | (ancestor[k] & child_r[k]); 
                                                  // Add tokens common to parent and either child to this node
    } else {                                     
      this[k] = this[k] | ancestor[k];            // Add parent's tokens to this node's tokens
      if (child_q[k] == *inapp || child_r[k] == *inapp) {
        continue;
      } else {
        if (this[k] & *inapp) this[k] -= *inapp;  // Remove {-} from this node's tokens
        (pars[k])++;                              // Increase parsimony score by one
        (*w) += weight[k];                        // Increase parsimony score by one
      }
    }
  }
}

void fitch_uppass_ia(int *state, int *parent_of, int *children_of, int *n_rows, int *pars, int *n_node, double *weight, double *pvec, int *inapp) {
  int i;
  fitch_uproot_ia(&state[(parent_of[0]-1L) * (*n_rows)], // this_start, will become this_finish
               &state[(children_of[0      ]-1L) * (*n_rows)], // child q
               &state[(children_of[*n_node]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0]-1L]);
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (i = 0; i < (*n_node)-1L; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // state and state are stored as [0 * n_rows] = states[,1], [1 * n_rows] = states[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    fitch_upnode_ia(
    //  The position of node 12 in the state array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &state[(parent_of[0] /*+1L -1L*/ + i) * (*n_rows)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need state[11-1]
    &state[(parent_of[i]-1L) * (*n_rows)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &state[(children_of[i + 1L]-1L) * (*n_rows)], // child q
    &state[(children_of[i + 1L + *n_node]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0] /*+1L -1L*/ + i]);
  }
}

SEXP FITCHUPIA(SEXP dat, SEXP n_transform_series, SEXP parent_of, SEXP children_of, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {
  int *state, *n_rows=INTEGER(n_transform_series), m=INTEGER(max_node)[0], i;   
  double *pvtmp;
  SEXP DATA, pars, pvec, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 3L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pvec = allocVector(REALSXP, m));
  PROTECT(DATA = allocMatrix(INTSXP, n_rows[0], m));
  pvtmp = REAL(pvec);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  state = INTEGER(DATA);
  for(i=0; i<(*n_rows * m); i++) state[i] = INTEGER(dat)[i];
  
  fitch_uppass_ia(state, INTEGER(parent_of), INTEGER(children_of), n_rows, INTEGER(pars), INTEGER(n_node), REAL(weight), pvtmp, INTEGER(inapp));
  
  SET_VECTOR_ELT(RESULT, 0, pars);
  SET_VECTOR_ELT(RESULT, 1, DATA);
  SET_VECTOR_ELT(RESULT, 2, pvec);
  UNPROTECT(4);
  return(RESULT); 
}


void fitch_downnode(int *dat1, int *dat2, int *n_rows, int *pars, double *weight, int *inapp, double *w, int *need_uppass) {
  int k, tmp;
  for (k = 0; k < (*n_rows); k++) {
    tmp = dat1[k] & dat2[k];
    if (tmp) {// Tokens in Common
      if ((tmp == (*inapp))) {  // Is - the only common token?
        need_uppass[k] = 1L;
        if ((dat1[k] | dat2[k]) != (*inapp)) { // At least one child has an applicable token        
          if ((dat1[k] != (*inapp)) && (dat2[k] != (*inapp))) { // ... Do both children have an applicable token?
            tmp = dat1[k] | dat2[k];
          }
        }
      }
    } else {
      tmp = dat1[k] | dat2[k];
      if ((dat1[k] == *inapp) || (dat2[k] == *inapp)) { // One child's only possible token {-}
      } else {
        if (tmp & (*inapp)) { // Only delete inapplicable token if it is present
          tmp = tmp - (*inapp);
        }
        (pars[k])++;
        (*w) += weight[k];
      }
    }
    dat1[k] = tmp;
  }
}

void fitch_downpass(int *dat, int *n_rows, int *pars, int *parent, int *child, int *n_edge, double *weight, int *inapp, double *pvec, double *pscore, int *need_uppass) {
  int i, ni, k;
  ni = 0;
  for (i = 0; i < *n_edge; i++) {
    if (ni == parent[i]){
      pvec[ni-1L] += pvec[child[i]-1L];
      fitch_downnode(&dat[(ni-1L) * (*n_rows)], &dat[(child[i]-1L) * (*n_rows)], n_rows, pars, weight, inapp, &pvec[(ni-1L)], need_uppass);
    } else {
      ni = parent[i];
      pvec[ni-1L] += pvec[child[i]-1L];        
      for (k = 0; k < (*n_rows); k++) dat[(ni-1L)*(*n_rows) + k] = dat[(child[i]-1L)*(*n_rows) + k];                     
    }
  }
  pscore[0] = pvec[ni-1];
}

SEXP FITCHDOWN(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP n_edge, SEXP weight, SEXP mx, SEXP q, SEXP inapp) {   
  int *data, *n_rows=INTEGER(nrx), *need_uppass, m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
  double *pvtmp;
  SEXP RESULT, pars, pscore, DAT, pvec, need_up;
  PROTECT(RESULT = allocVector(VECSXP, 5L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  PROTECT(need_up = allocVector(INTSXP, *n_rows));
  pvtmp = REAL(pvec);
  need_uppass = INTEGER(need_up);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) {INTEGER(pars)[i] = 0L; need_uppass[i] = 0L;}
  REAL(pscore)[0]=0.0;
  data = INTEGER(DAT);
  for(i=0; i<(*n_rows * n); i++) data[i] = INTEGER(dat)[i];
  
  fitch_downpass(data, n_rows, INTEGER(pars), INTEGER(parent), INTEGER(child), INTEGER(n_edge), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore), need_uppass);
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  SET_VECTOR_ELT(RESULT, 4, need_up);
  UNPROTECT(6);
  return(RESULT);
}

void fitch_uproot(int *this, int *child_q, int *child_r, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, ancestor_k;
  for (k = 0; k < (*n_rows); k++) { // Next TS
    ancestor_k = this[k];
    if ((ancestor_k & (*inapp)) && (ancestor_k != (*inapp))) ancestor_k -= (*inapp);  // Remove {-} from parent node: Hennig's Auxiliary Principle
    if (this[k] & (*inapp)) {       // Node has an inapplicable token?
      if (ancestor_k == *inapp) {   // Parent is {-}
        this[k] = *inapp;           // Set node's tokens to {-}
        continue;
      } 
      else 
      if (this[k] != (*inapp)) {    // Node has an applicable token?
        this[k] -= (*inapp);        // Remove {-} from this node's tokens
        if (child_q[k] == (*inapp) || child_r[k] == (*inapp)                      // One child's only possible token is {-}
        || ((child_q[k] & child_r[k]) && ((child_q[k] & child_r[k]) != (*inapp))) // Children have tokens in common, excluding {-}
        ) {} else {
          (pars[k])++; (*w) += weight[k];  // Increase parsimony score by one
          this[k] = this[k] | (ancestor_k & child_r[k]) | (ancestor_k & child_q[k]);
          continue;                 // Next node
        }
      } else if (ancestor_k & child_r[k] & child_q[k] & *inapp) {
        continue;                   // Parent and both children have inapplicable token
      } else if (ancestor_k != (*inapp) && (child_q[k] != (*inapp) || child_r[k] != (*inapp))) { // Parent + one child have applicable token
        if (!((((child_q[k] | child_r[k]) & ancestor_k) | (*inapp)) != (*inapp))) { // Parent has no applicable token in common with either child
          (pars[k])++; (*w) += weight[k];  // Increase parsimony score by one
          this[k] = child_q[k] | child_r[k] | ancestor_k; // Set this node's tokens to the tokens present in the parent or either child...
        } else {
          this[k] = (child_q[k] | child_r[k]) & ancestor_k; // Set this node's tokens to those in common...
        }
        if (this[k] & *inapp) this[k] -= *inapp;             // Excluding inapplicable
        continue;                   // Next node
      }
    }
    if ((ancestor_k & this[k]) == ancestor_k) { // All parent node's tokens among this node's possible tokens
      this[k] = ancestor_k;                      // Set this node's tokens to parent's tokens
    } else if (child_q[k] & child_r[k]) {         // Children have tokens in common
      this[k] = this[k] | (ancestor_k & child_q[k]) | (ancestor_k & child_r[k]); 
                                                  // Add tokens common to parent and either child to this node
    } else {                                     
      this[k] = this[k] | ancestor_k;            // Add parent's tokens to this node's tokens
    }
  }
}

void fitch_upnode(int *this, int *ancestor, int *child_q, int *child_r, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k;
  for (k = 0; k < (*n_rows); k++) { // Next node in preorder
    if (this[k] & (*inapp)) {       // Node has an inapplicable token?
      if (ancestor[k] == *inapp) {  // Parent is {-}
        this[k] = *inapp;           // Set node's tokens to {-}
        continue;
      } 
      else if (this[k] != (*inapp)) {    // Node has an applicable token?
        this[k] -= (*inapp);        // Remove {-} from this node's tokens
        if (child_q[k] == (*inapp) || child_r[k] == (*inapp)                      // One child's only possible token is {-}
        || ((child_q[k] & child_r[k]) && ((child_q[k] & child_r[k]) != (*inapp))) // Children have tokens in common, excluding {-}
        ) {} else {
          (pars[k])++; (*w) += weight[k];  // Increase parsimony score by one
          this[k] = this[k] | (ancestor[k] & child_r[k]) | (ancestor[k] & child_q[k]);
          continue;                 // Next node
        }
      } else if (ancestor[k] & child_r[k] & child_q[k] & *inapp) {
        continue;                   // Parent and both children have inapplicable token
      } else if (ancestor[k] != (*inapp) && (child_q[k] != (*inapp) || child_r[k] != (*inapp))) { // Parent + one child have applicable token
        if (!((((child_q[k] | child_r[k]) & ancestor[k]) | (*inapp)) != (*inapp))) { // Parent has no applicable token in common with either child
          (pars[k])++; (*w) += weight[k];  // Increase parsimony score by one
          this[k] = child_q[k] | child_r[k] | ancestor[k]; // Set this node's tokens to the tokens present in the parent or either child...
        } else {
          this[k] = (child_q[k] | child_r[k]) & ancestor[k]; // Set this node's tokens to those in common...
        }
        if (this[k] & *inapp) this[k] -= *inapp;             // Excluding inapplicable
        continue;                   // Next node
      }
    }
    if ((ancestor[k] & this[k]) == ancestor[k]) { // All parent node's tokens among this node's possible tokens
      this[k] = ancestor[k];                      // Set this node's tokens to parent's tokens
    } else if (child_q[k] & child_r[k]) {         // Children have tokens in common
      this[k] = this[k] | (ancestor[k] & child_q[k]) | (ancestor[k] & child_r[k]); 
                                                  // Add tokens common to parent and either child to this node
    } else {                                     
      this[k] = this[k] | ancestor[k];            // Add parent's tokens to this node's tokens
    }
    if (this[k] & *inapp) {
      continue;
    } else {
      if ((child_q[k] & *inapp) && child_q[k] != *inapp) child_q[k] -= *inapp;
      if ((child_r[k] & *inapp) && child_r[k] != *inapp) child_r[k] -= *inapp;
    }
  }
}

void fitch_uppass(int *state, int *parent_of, int *children_of, int *n_rows, int *pars, int *n_node, double *weight, double *pvec, int *inapp, double *pscore) {
  int i;
  fitch_uproot(&state[(parent_of[0]-1L) * (*n_rows)], // this_start, will become this_finish
               &state[(children_of[0      ]-1L) * (*n_rows)], // child q
               &state[(children_of[*n_node]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0]-1L]);
  pscore[0] += pvec[parent_of[0] - 1L];
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (i = 0; i < (*n_node)-1L; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // state and state are stored as [0 * n_rows] = states[,1], [1 * n_rows] = states[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    fitch_upnode(
    //  The position of node 12 in the state array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &state[(parent_of[0] /*+1L -1L*/ + i) * (*n_rows)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need state[11-1]
    &state[(parent_of[i]-1L) * (*n_rows)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &state[(children_of[i + 1L]-1L) * (*n_rows)], // child q
    &state[(children_of[i + 1L + *n_node]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0] /*+1L -1L*/ + i]);
    pscore[0] += pvec[parent_of[0] /*+1L -1L*/ + i];
  }
}

SEXP FITCHUP(SEXP dat, SEXP n_transform_series, SEXP parent_of, SEXP children_of, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {
  int *state, *n_rows=INTEGER(n_transform_series), m=INTEGER(max_node)[0], i;   
  double *pvtmp;
  SEXP DATA, pars, pscore, pvec, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(pvec = allocVector(REALSXP, m));
  PROTECT(DATA = allocMatrix(INTSXP, n_rows[0], m));
  pvtmp = REAL(pvec);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0] = 0.0;
  state = INTEGER(DATA);
  for(i=0; i<(*n_rows * m); i++) state[i] = INTEGER(dat)[i];
  
  fitch_uppass(state, INTEGER(parent_of), INTEGER(children_of), n_rows, INTEGER(pars), INTEGER(n_node), REAL(weight), pvtmp, INTEGER(inapp), REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DATA);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}