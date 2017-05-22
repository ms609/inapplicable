#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void app_fitch_downnode(int *this, int *left, int *right, int start_char, int *end_char, int *pars, double *weight, double *w) {
  int i;
  for (i = start_char; i < (*end_char); i++) {
    if (left[i] & right[i]) {
      this[i] = left[i] & right[i];
    }
    else {
      this[i] = left[i] | right[i];
      (pars[i])++; (*w) += weight[i]; // Add one to tree length
    }
  }
}

void app_fitch_downpass
(int *dat, int start_char, int *n_char, int *pars, int *parent, int *child, int *n_edge, double *weight, int *inapp, double *pvec, double *pscore) {
  int parent_i = 0;
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    pvec[parent_i -1] += pvec[child[i]-1] + pvec[child[i+1]-1];
    fitch_app_downnode(&dat[(parent_i-1)* (*n_char)], &dat[(child[i]-1) * (*n_char)], &dat[(child[i+1]-1) * (*n_char)], start_char, n_char, pars, weight, &pvec[(parent_i-1)]);
  }
  pscore[0] = pvec[parent_i-1];
}

void app_upnode
(int *this, int *app, int *ancestor, int *left, int *right, int inapp, int n_char) {
  int i;
  for (i = 0; i < n_char; i++) {
    if (this[i] & inapp) {
      if (this[i] & ~inapp) {
        if (anc[i] == inapp) {
          app[i] = inapp;
        }
        else {
          app[i] = this[i] & ~inapp;
        }
      }
      else {
        if (anc[i] == inapp) {
          app[i] = inapp;
        }
        else {
          if ((left[i] | right[i]) & ~inapp) {
            app[i] = ((left[i] | right[i]) & ~inapp);
          }
          else {
            app[i] = inapp;
          }
        }
      }
    }
    else {
      app[i] = this[i];
    }      
    assert(app[i]);
  }
}

void app_uproot
(int *this, int inapp, int n_char) {
  for (int i = 0; i < n_char; i++) this[i] = ((this[i] == inapp) ? inapp : this[i] & ~inapp); // Assume applicable if ambiguous.
}

void inapp_first_uppass
(int *app, int n_char, int *parent_of, int *children_of, int *n_node) {
  app_uproot(&app[(parent_of[0]-1) * n_char], n_char);
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
    &app[(children_of[i + 1] -1) * (*n_char)], // left child
    &app[(children_of[i + 1 + *n_node] -1) * (*n_char)], // right child
    n_char);
  }
}

void inapp_first_downnode
(int *dat, int *left, int *right, int *this_acts, int *l_acts, int *r_acts, int inapp, int n_char) {
  int i, tmp;
  for (i = 0; i < nchars; ++i) {
    if ((temp = (left[i] & right[i]))) {
      this[i] = temp;
      
      if (temp == inapp) {
        if ((left[i] & ~inapp) && (right[i] & ~inapp)) {
          this[i] = (left[i] | right[i]);
        }
      }
    }
    else {
      this[i] = (left[i] | right[i]);
      
      if ((left[i] & ~inapp) && (right[i] & ~inapp)) {
          this[i] = this[i] & ~inapp;
      }
    }
    
    this_acts[i] = (l_acts[i] | r_acts[i]) & ~inapp;
    
    assert(this[i]);
  }
}

void inapp_first_downpass
(int *app, int n_char, int *parent, int *child, int *n_edge) {
  for (int i = 0; i < *n_edge; i+=2) {
    inapp_first_downnode(&app[(parent[i]  - 1) * n_char],
                         &app[( child[i]  - 1) * n_char],
                         &app[( child[i+1]- 1) * n_char], n_char);
  }
}

void inapp_second_downnode
(int *this, int *this_app, int *left, int *right, int *this_acts, int *l_acts, int *r_acts
 int end_char, int inapp, int *pars, double *weight, double *w) {
  int temp;
  for (i = 0; i < end_char; ++i) {
    if (this_app[i] != inapp) {
      if ((temp = (left[i] & right[i]))) {
        if (temp & ~inapp) {
          this[i] = temp & ~inapp;
        } else {
          this[i] = temp;
        }
      }
      else {
        this[i] = (left[i] | right[i]) & ~inapp;
        
        if ((left[i] & ~inapp && right[i] & ~inapp)
        ||  (l_acts[i] && r_acts[i])) {
          (pars[i])++; (*w) += weight[i]; // Add one to tree length
        }
      }
    }
    
    this_acts[i] = (l_acts[i] | r_acts[i]) & ~inapp;

    assert(this[i]);
  }
  
  return steps;
}

void inapp_second_downpass
(int *dat, int *app, int *act, int *n_char, int *pars, int *parent, int *child, int *n_edge, double *weight, int *inapp, double *pvec, double *pscore) {
  int parent_i = 0;
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    pvec[parent_i -1] += pvec[child[i]-1] + pvec[child[i+1]-1];
    fitch_app_downnode(
      &dat[(parent_i-1)* (*n_char)],
      &app[(parent_i-1)* (*n_char)],
      &dat[(child[i]-1) * (*n_char)], 
      &dat[(child[i+1]-1) * (*n_char)], 
      &act[(parent_i-1)* (*n_char)],
      &act[(child[i]-1) * (*n_char)], 
      &act[(child[i+1]-1) * (*n_char)], 
      *n_char, *inapp, pars, weight, &pvec[(parent_i-1)]
    );
  }
  pscore[0] = pvec[parent_i-1];  
}

void inapp_second_upnode
(int *this, int *left, int *right, int *this_acts, int *l_acts, int *r_acts
 int end_char, int inapp, int *pars, double *weight, double *w) {
  for (i = 0; i < nchars; ++i) {    
    if (this[i] & ~inapp) {
      if (anc[i] & ~inapp) {
        if ((anc[i] & this[i]) == anc[i]) {
          this[i] = anc[i] & this[i];
        } else {
          if (left[i] & right[i]) {
            this[i] = (this[i] | (anc[i] & left[i] & right[i]));
          }
          else {
            if ((left[i] | right[i]) & inapp) {
              if ((left[i] | right[i]) & anc[i]) {
                this[i] = anc[i];
              } else {
                this[i] = (left[i] | right[i] | anc[i]) & inapp;
              }
            } else {
              this[i] = this[i] | anc[i];
              if ((anc[i] & this[i]) == anc[i]) {
                this[i] = anc[i] & this[i];
              }
            }
          }
        }
      }
    }
    else {
      if (l_acts[i] && r_acts[i]) {
        (pars[i])++; (*w) += weight[i]; // Add one to tree length
      }
    }
    assert(this[i]);
  }

}

void inapp_second_uppass
(int *dat, int *act, int start_char, int *n_char, int *pars, int *parent, int *child, int *n_edge, double *weight, int *inapp, double *pvec, double *pscore) {
  app_uproot(&app[(parent_of[0]-1) * n_char], n_char);
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
    &app[(children_of[i + 1] -1) * (*n_char)], // left child
    &app[(children_of[i + 1 + *n_node] -1) * (*n_char)], // right child
    n_char);
  }
  
}

SEXP MORPHYFITCH(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP parent_of, SEXP children_of, SEXP n_edge, SEXP n_node, SEXP weight, SEXP max_node, SEXP n_tip, SEXP inapp, SEXP inapp_chars) {   
  int *data, *actives, *n_char=INTEGER(nrx), *inappl=INTEGER(inapp), *appl,  m=INTEGER(max_node)[0], 
      i, n=INTEGER(n_tip)[0], first_applicable=INTEGER(inapp_chars)[0];
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
    actives[i] = 0;
  }
  
  inapp_first_downpass(appl, first_applicable, INTEGER(parent), INTEGER(child), INTEGER(n_edge));
  inapp_first_uppass  (appl, first_applicable, INTEGER(parent), INTEGER(child), INTEGER(n_edge));
  // TODO!: inapp_update_tips   (appl, first_applicable, INTEGER(parent), INTEGER(child), INTEGER(n_edge)); 
  inapp_second_downpass(data, appl, actives, n_char, INTEGER(pars), INTEGER(parent), INTEGER(child), INTEGER(n_edge), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  inapp_second_uppass  (data, actives, n_char, INTEGER(pars), INTEGER(parent), INTEGER(child), INTEGER(n_edge), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  
  app_fitch_downpass(data, first_applicable, n_char, INTEGER(pars), INTEGER(parent), INTEGER(child), INTEGER(n_edge), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore)); // No need for an up-pass: all scoring on way down.

  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  SET_VECTOR_ELT(RESULT, 4, APPL);
  UNPROTECT(6);
  return(RESULT);
}