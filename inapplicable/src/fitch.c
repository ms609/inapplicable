#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void fitchdown(int *dat1, int *dat2, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, tmp, noin1, noin2;
  for (k = 0; k < (*n_rows); k++) {
    tmp = dat1[k] & dat2[k];
    if (tmp // Tokens in Common
      && ((*inapp) & tmp) // && Inapplicable in common 
      && (noin1 = dat1[k] - (*inapp)) && (noin2 = dat2[k] - (*inapp)) // && Apart from inapplicable, tokens
      && !((noin1) & (dat2[k] - (*inapp))) // && Apart from inapplicable, no tokens in common
    ){
      pars[k] = -1; // Flag that we'll need to perform an uppass on transformation series k
      dat1[k] = dat1[k] | dat2[k];
      continue;
    } else if (!tmp){
      tmp = dat1[k] | dat2[k];
      if ((dat1[k] == *inapp) || (dat2[k] == *inapp)) { // Do not increment pscore
      } else {
        if ((*inapp) & tmp) { // Only delete inapplicable token if it is present
          tmp = tmp - (*inapp);
        }
        (pars[k])++;
        (*w) += weight[k];  // TODO: does this need to be in a temporary w_tmp, only added to *w when we're sure we won't be fed a 'continue'?
      }
    }
    dat1[k] = tmp;
  }
}

void fitchbridge(int *dat, int *n_rows, int *pars, int *node, int *edge, int *nl, double *weight, int *inapp, double *pvec, double *pscore) {
  int i, ni, k;
  ni = 0;
  for (i=0; i< *nl; i++) {
    if (ni == node[i]){
      pvec[ni-1L] += pvec[edge[i]-1L];
      fitchdown(&dat[(ni-1L) * (*n_rows)], &dat[(edge[i]-1L) * (*n_rows)], n_rows, pars, weight, inapp, &pvec[(ni-1L)]);
    } else {
      ni = node[i];
      pvec[(ni-1L)] += pvec[(edge[i]-1L)];         
      for(k = 0; k < (*n_rows); k++) dat[(ni-1L)*(*n_rows) + k] = dat[(edge[i]-1L)*(*n_rows) + k];                     
    }
  }
  pscore[0]=pvec[ni-1];
}

SEXP FITCHI(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q, SEXP inapp) {   
  int *data, *n_rows=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
  double *pvtmp;
  SEXP DAT, pars, pvec, pscore, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  pvtmp = REAL(pvec);
  data = INTEGER(DAT);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0]=0.0;
  for(i=0; i<(*n_rows * n); i++) data[i] = INTEGER(dat)[i];
  
  fitchbridge(data, n_rows, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}

void fitch_upnode(int *this, int *ancestor, int *childq, int *childr, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, final;
  for (k = 0; k < (*n_rows); k++) { // Next TS
    if (!(ancestor[k] & (*inapp))   // Does parent lack {-}?
      && (this[k] & (*inapp))       // Does this node have {-}?
      && (this[k] - (*inapp))       // Does this node have an applicable token?
    ) {
      final -= (*inapp);            // Remove {-} from this node's tokens
      if (childq[k] == (*inapp) || childr[k] == (*inapp)                  // One child's only possible token is {-}
      || ((childq[k] & childr[k]) && ((childq[k] & childr[k]) != (*inapp))) // Children have tokens in common, excluding {-}
      ) {} else {
        pars[k] += 1;
        this[k] = (ancestor[k] & childr[k]) | (ancestor[k] & childq[k]) | final;
        continue;
      }
    }
    if ((ancestor[k] & final) == ancestor[k]) { // All parent node's tokens among this node's possible tokens
      this[k] = ancestor[k];                  // Set this node's tokens to parent's tokens
    } else if (childq[k] & childr[k]) {       // Children have tokens in common
      this[k] = final | ancestor[k];          // Add parent's tokens to this node's tokens
    } else {                                  // Add tokens common to parent and either child to this node
      this[k] = (ancestor[k] & childr[k]) | (ancestor[k] & childq[k]) | final; 
    }
  }
}

void fitch_uppass(int *state, int *n_rows, int *pars, int *parent_of, int *child_of, int *n_node, double *weight, int *inapp, double *pvec, double *pscore) {
  int i;
  fitch_upnode(&state[(parent_of[0]) * (*n_rows)], // this_start, will become this_finish
  &state[(parent_of[0]) * (*n_rows)], // ancestor_finish
  &state[(child_of[0]-1L) * (*n_rows)], &state[(child_of[1L]-1L) * (*n_rows)], // childq_start, childr_start
  n_rows, pars, weight, inapp, &pvec[0]); // #TODO! checkL pvec 0 or pvec[root]?
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (i = 0; i < (*n_node)-1L; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // child_of  is stored as 0L, 0R, 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR
    // state and state are stored as [0 * n_rows] = states[,1], [1 * n_rows] = states[,2], ....
    // Worked examples assume that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    fitch_upnode(
    //  The position of node 12 in the state and state array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &state[(parent_of[0] + i) * (*n_rows)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need state[11-1]
    &state[(parent_of[i]-1L) * (*n_rows)], // ancestor_finish
    //  To find the number of node 12's children we look in child_of[node12.index]
    //    child_of[0, 1] are the two children of node [root + i] = 12
    //    node12.index = i = 0; child_of[0*2] = Q; child_of[0*2 + 1] = R
    &state[(child_of[i * 2L]-1L) * (*n_rows)], &state[(child_of[(i * 2L) + 1L]-1L) * (*n_rows)], // childq_start, childr_start
    n_rows, pars, weight, inapp, &pvec[i]);
    // # TODO : Worry about pvec
  }
  // pscore[0]=pvec[ni-1];
}

SEXP FITCHUP(SEXP dat, SEXP n_transform_series, SEXP parent_of, SEXP child_of, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {
  int *state, *n_rows=INTEGER(n_transform_series), m=INTEGER(max_node)[0], i, n=INTEGER(n_tip)[0];   
  double *pvtmp;
  SEXP DAT, pars, pvec, pscore, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  pvtmp = REAL(pvec);
  state = INTEGER(DAT);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0]=0.0;
  for(i=0; i<(*n_rows * n); i++) state[i] = INTEGER(dat)[i];
  state = state;
  
  fitch_uppass(state, n_rows, INTEGER(pars), INTEGER(parent_of), INTEGER(child_of), INTEGER(n_node), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}