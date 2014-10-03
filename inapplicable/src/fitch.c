#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void fitch_app_downnode(int *app, int *this, int *child_q, int *child_r, int *n_char, int *pars, double *weight, int *inapp, double *w) {
  int tmp;
  for (int k = 0; k < (*n_char); k++) {
    if (app[k] == -1) {
      this[k] = *inapp;
    } else if (child_q[k] == *inapp) {
      this[k] = child_r[k] & ~*inapp;
    } else if (child_r[k] == *inapp) {
      this[k] = child_q[k] & ~*inapp;
    } else if ((tmp = (child_q[k] & child_r[k] & ~*inapp))) { // Applicable tokens in Common
      this[k] = tmp;
    } else {
      this[k] = (child_q[k] | child_r[k]) & ~*inapp;
      (pars[k])++; (*w) += weight[k]; // Add one to tree length
    }
  }
}

void fitch_app_downpass(int *dat, int *app, int *n_char, int *pars, int *parent, int *child, int *n_edge, double *weight, int *inapp, double *pvec, double *pscore) {
  int parent_i = 0;
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    pvec[parent_i -1] += pvec[child[i] -1] + pvec[child[i+1] -1];
    fitch_app_downnode(&app[(parent_i -1)*(*n_char)], &dat[(parent_i -1)* (*n_char)], &dat[(child[i]-1) * (*n_char)], &dat[(child[i+1]-1) * (*n_char)], n_char, pars, weight, inapp, &pvec[(parent_i -1)]);
  }
  pscore[0] = pvec[parent_i-1];
}

void app_upnode(int *this, int *ancestor, int *child_q, int *child_r, int *n_char) {
  for (int k = 0; k < (*n_char); k++) this[k] = ((child_q[k] + child_r[k] + ancestor[k] >= 0) ? 1 : -1);
}

void app_uproot(int *this, int *n_char) {
  for (int k = 0; k < (*n_char); k++) this[k] = ((this[k] > -1) ? 1 : -1); // Assume applicable if ambiguous.
}

void app_uppass(int *app, int *n_char, int *parent_of, int *children_of, int *n_node) {
  app_uproot(&app[(parent_of[0]-1) * (*n_char)], n_char);
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (int i = 0; i < (*n_node)-1; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // app and app are stored as [0 * n_char] = apps[,1], [1 * n_char] = apps[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    app_upnode(
    //  The position of node 12 in the app array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &app[(parent_of[0] + i + 1 -1) * (*n_char)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need app[11-1]
    &app[(parent_of[i]-1) * (*n_char)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &app[(children_of[i + 1]-1) * (*n_char)], // child q
    &app[(children_of[i + 1 + *n_node]-1) * (*n_char)], // child r
    n_char);
  }
}

void app_downnode(int *dat, int *child_q, int *child_r, int *n_char) {
  int tmp;
  for (int k = 0; k < (*n_char); k++) {
    dat[k] = ((tmp = (child_q[k] + child_r[k])) ? (tmp > 0 ? 1 : -1) : 0);
  }
}

void app_downpass(int *app, int *n_char, int *parent, int *child, int *n_edge) {
  for (int i = 0; i < *n_edge; i+=2) {
    app_downnode(&app[(parent[i]  - 1) * (*n_char)],
                 &app[( child[i]  - 1) * (*n_char)],
                 &app[( child[i+1]- 1) * (*n_char)], n_char);
  }
}

SEXP FITCHINAPP(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP parent_of, SEXP children_of, SEXP n_edge, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp) {   
  int *data, *n_char=INTEGER(nrx), *inappl=INTEGER(inapp), *appl,  m=INTEGER(max_node)[0], i, n=INTEGER(n_tip)[0];
  double *pvtmp;
  SEXP RESULT, pars, pscore, DAT, pvec, APPL;
  PROTECT(RESULT = allocVector(VECSXP, 5));
  PROTECT(pars = allocVector(INTSXP, *n_char));
  PROTECT(pscore = allocVector(REALSXP, 1));
  PROTECT(DAT = allocMatrix(INTSXP, n_char[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  PROTECT(APPL = allocMatrix(INTSXP, n_char[0], m));
  pvtmp = REAL(pvec);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_char; i++) {INTEGER(pars)[i] = 0;}
  REAL(pscore)[0]=0.0;
  data = INTEGER(DAT);
  appl = INTEGER(APPL);
  for(i=0; i<(*n_char * n); i++) {
    data[i] = INTEGER(dat)[i];
    appl[i] = ((data[i] & ~*inappl) ? 1 : 0) - ((data[i] & *inappl) ? 1 : 0);
  }
  
  app_downpass(appl, n_char, INTEGER(parent), INTEGER(child), INTEGER(n_edge));
  app_uppass(  appl, n_char, INTEGER(parent_of), INTEGER(children_of), INTEGER(n_node));
  fitch_app_downpass(data, appl, n_char, INTEGER(pars), INTEGER(parent), INTEGER(child), INTEGER(n_edge), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  SET_VECTOR_ELT(RESULT, 4, APPL);
  UNPROTECT(6);
  return(RESULT);
}