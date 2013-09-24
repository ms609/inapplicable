#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void fitch_downnode(int *dat1, int *dat2, int *n_rows, int *pars, double *weight, int *inapp, double *w, int *need_uppass) {
  int k, tmp, noin1, noin2;
  for (k = 0; k < (*n_rows); k++) {
    tmp = dat1[k] & dat2[k];
    if (tmp // Tokens in Common
      && ((*inapp) & tmp) // Inapplicable in common 
      && (noin1 = dat1[k] - (*inapp)) && (noin2 = dat2[k] - (*inapp)) // Apart from inapplicable, tokens
      && !((noin1) & (dat2[k] - (*inapp))) // Apart from inapplicable, no tokens in common
    ) {
      need_uppass[k] = 1L;
      dat1[k] = dat1[k] | dat2[k];
    } else if (!tmp){
      tmp = dat1[k] | dat2[k];
      if ((dat1[k] == *inapp) || (dat2[k] == *inapp)) { // Do not increment pscore
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

SEXP FITCHI(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP n_edge, SEXP weight, SEXP mx, SEXP q, SEXP inapp) {   
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
  data = INTEGER(DAT);
  need_uppass = INTEGER(need_up);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) {INTEGER(pars)[i] = 0L; need_uppass[i] = 0L;}
  REAL(pscore)[0]=0.0;
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

void fitch_upnode(int *this, int *ancestor, int *childq, int *childr, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, final;
  for (k = 0; k < (*n_rows); k++) { // Next TS
    final = this[k];
    if (!(ancestor[k] & (*inapp))   // Does parent lack {-}?
      && (this[k] & (*inapp))       // Does this node have {-}?
      && (this[k] - (*inapp))       // Does this node have an applicable token?
    ) {
      final -= (*inapp);            // Remove {-} from this node's tokens
      if (childq[k] == (*inapp) || childr[k] == (*inapp)                  // One child's only possible token is {-}
      || ((childq[k] & childr[k]) && ((childq[k] & childr[k]) != (*inapp))) // Children have tokens in common, excluding {-}
      ) {} else {
        (pars[k])++;
        (*w) += weight[k];
        this[k] = (ancestor[k] & childr[k]) | (ancestor[k] & childq[k]) | final;
        *w = 1;
        continue;
      }
    }
    if ((ancestor[k] & final) == ancestor[k]) { // All parent node's tokens among this node's possible tokens
      *w = ancestor[k] * 1000;
      this[k] = ancestor[k];                  // Set this node's tokens to parent's tokens
    } else if (childq[k] & childr[k]) {       // Children have tokens in common
      *w = 3;
      this[k] = final | ancestor[k];          // Add parent's tokens to this node's tokens
    } else {                                  // Add tokens common to parent and either child to this node
      *w = 5;
      this[k] = (ancestor[k] & childq[k]) | (ancestor[k] & childr[k]) | final; 
    }
  }
}

void fitch_uppass(int *state, int *parent_of, int *children_of, int *n_rows, int *pars, int *n_node, double *weight, double *pvec, int *inapp, double *pscore) {
  int i;
  fitch_upnode(&state[(parent_of[0]-1L) * (*n_rows)], // this_start, will become this_finish
               &state[(parent_of[0]-1L) * (*n_rows)], // ancestor = this_start
               &state[(children_of[0 ]-1L) * (*n_rows)], // child q
               &state[(children_of[1L]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0]-1L]);
  pscore[0] += pvec[parent_of[0]];
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (i = 0; i < (*n_node)-1L; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of  is stored as 0L, 0R, 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR
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
    //    children_of[0, 1] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &state[(children_of[(i * 2L) + 2L]-1L) * (*n_rows)], // child q
    &state[(children_of[(i * 2L) + 3L]-1L) * (*n_rows)], // child r
    n_rows, pars, weight, inapp, &pvec[parent_of[0] /*+1L -1L*/ + i]);
    pscore[0] += pvec[parent_of[0] /*+1L -1L*/ + i];
  }
}

SEXP FITCHUP(SEXP dat, SEXP n_transform_series, SEXP parent_of, SEXP children_of, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {
  int *state, *n_rows=INTEGER(n_transform_series), m=INTEGER(max_node)[0], i, n=INTEGER(n_tip)[0];   
  double *pvtmp;
  SEXP DAT, pars, pscore, pvec, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(pvec = allocVector(REALSXP, m));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  state = INTEGER(DAT);
  pvtmp = REAL(pvec);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0] = 0.0;
  for(i=0; i<(*n_rows * n); i++) state[i] = INTEGER(dat)[i];
  
  fitch_uppass(state, INTEGER(parent_of), INTEGER(children_of), n_rows, INTEGER(pars), INTEGER(n_node), REAL(weight), pvtmp, INTEGER(inapp), REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}