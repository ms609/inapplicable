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
        (*w) += weight[k];
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

void fitch_uppass(int *this_start, int *this_finish, int *ancestor_finish, int *childq_start, int *childr_start, int *n_rows, int *pars, double *weight, int *inapp, double *w) {
  int k, tmp_finish;
  for (k = 0; k < (*n_rows); k++) {
    tmp_finish = dat1[k] & dat2[k]; // Intersect
    if (tmp_finish != ancestor_finish[k]) {
      if (childq_start & childr_start) {
        tmp_finish <- ((childq_start | childr_start) & ancestor_finish) | this_start
      } else {
        tmp_finish <- this_start | ancestor_finish // Union
      }
    }
    this_finish[k] = tmp_finish;
}

void fitchupbridge(int *start_state, int *n_rows, int *pars, int *parent_of, int *child_of, int *n_node, double *weight, int *inapp, double *pvec, double *pscore) {
  int i, prev_parent, k;
  prev_parent = 0;
  &finsh_state[(parent_of[0]-1L * (*n_rows)] = &start_state[(parent_of[0]-1L * (*n_rows)]; //Node 0, i.e. root
  for (i = 1; i < *n_node; i++) { 
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // child_of  is stored as 0L, 0R, 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR
    // start_state and finish_state are stored as [0 * n_rows] = states[,1], [1 * n_rows] = states[,2], ....
    fitch_uppass(&start_state[(*n_node + i + 2L) * (*n_rows)], &finish_state[(*n_node + i + 2L) * (*n_rows)], &start_state[(parent_of[i-1L]-1L) * (*n_rows)], &start_state[(child_of[2L*i-1L]-1L) * (*n_rows)], &start_state[(child_of[2L*i]-1L) * (*n_rows)], nr, pars, weight, inapp, &pvec[(prev_parent-1L)]);
    /*
    if (prev_parent == parent[i]) {
      // pvec[prev_parent-1L] += pvec[child[i]-1L]; // work out what this is doing later
      fitchuppass(&dat[(prev_parent-1L) * (*n_rows)], &dat[(child[i]-1L) * (*n_rows)], nr, pars, weight, inapp, &pvec[(prev_parent-1L)]);
    } else {
      prev_parent = parent[i];
      // pvec[(prev_parent-1L)] += pvec[(child[i]-1L)]; // work out what this is doing later
      // for(k = 0; k < (*n_rows); k++) dat[(prev_parent-1L)*(*n_rows) + k] = dat[(child[i]-1L)*(*n_rows) + k]; // work out what this is doing later
    }
  */
  }
  pscore[0]=pvec[ni-1];
}

SEXP FITCHUP(SEXP dat, SEXP n_transform_series, SEXP parent_of, SEXP child_of, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {   
  int *data, *n_rows=INTEGER(n_transform_series), m=INTEGER(max_node)[0], i, n=INTEGER(n_tip)[0];   
  double *pvtmp;
  SEXP DAT, pars, pvec, pscore, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *n_rows));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  pvtmp = REAL(pvec);
  data = INTEGER(DAT);
  finish_state = INTEGER(DAT);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_rows; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0]=0.0;
  for(i=0; i<(*n_rows * n); i++) data[i] = INTEGER(dat)[i];
  
  fitchupbridge(data, finish_state, n_rows, INTEGER(pars), INTEGER(parent_of), INTEGER(child_of), INTEGER(n_node), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}