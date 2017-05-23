#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

// Note: Where possible pass pointers rather than their values to conserve stack space.

void app_fitch_downnode
(int *this, int *left, int *right, int *start_char, int *end_char, 
 int *pars, double *weight, double *w) {
  int i;
  for (i = *start_char; i < (*end_char); i++) {
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
(int *dat, int *parent, int *child, int *start_char, int *n_char, int *n_edge, 
 int *inapp, int *pars, double *weight, double *pvec) {
  int parent_i = 0;
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    pvec[parent_i -1] += pvec[child[i+1]-1] + pvec[child[i]-1];
    app_fitch_downnode(&dat[(parent_i  -1) * (*n_char)],
                       &dat[(child[i+1]-1) * (*n_char)], 
                       &dat[(child[i]  -1) * (*n_char)],
                       start_char, n_char, pars, weight, &pvec[(parent_i-1)]);
  }
}

void inapp_first_upnode
(int *this, int *app, int *ancestor, int *left, int *right,
 int *inapp, int *end_char) {
  int i;
  for (i = 0; i < *end_char; i++) {
    if (this[i] & (*inapp)) {
      if (this[i] & ~(*inapp)) {
        if (ancestor[i] == *inapp) {
          app[i] = *inapp;
        }
        else {
          app[i] = this[i] & ~(*inapp); // TODO posit: just ~(*inapp) will do
        }
      }
      else {
        if (ancestor[i] == *inapp) {
          app[i] = *inapp;
        }
        else {
          if ((left[i] | right[i]) & ~(*inapp)) {
            app[i] = ((left[i] | right[i]) & ~(*inapp)); // TODO posit: just ~(*inapp) will do
          }
          else {
            app[i] = *inapp;
          }
        }
      }
    }
    else {
      app[i] = this[i];
    }
  }
}

void inapp_first_root(int *this, int *app, int *inapp, int *end_char) {
  for (int i = 0; i < *end_char; i++) {
    if (this[i] != *inapp) {
      app[i] = this[i] & ~(*inapp); // TODO posit: ~(*inapp) will do?
    } else {
      app[i] = *inapp;
    }
  }
}

void inapp_first_uppass
(int *dat, int *app, int *parent_of, int *children_of, 
 int *end_char, int *n_char, int *n_node, int *inapp) {
  inapp_first_root(&dat[(parent_of[0]-1) * (*n_char)], 
                   &app[(parent_of[0]-1) * (*n_char)], inapp, end_char);
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (int i = 0; i < (*n_node) - 1; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // app and app are stored as [0 * n_char] = apps[,1], [1 * n_char] = apps[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    inapp_first_upnode(
      //  The position of node 12 in the dat / app array is:
      //    root.number (counting from 1) + i = i.e. position[11], the 12th position
      &dat[(parent_of[0] + i + 1 -1) * (*n_char)], // this_start
      &app[(parent_of[0] + i + 1 -1) * (*n_char)], // store applicability for future ref
      //  To find the number of node 12's parent we look in parent_of[node12.index]
      //    parent_of[0] is the parent of node [root + i] = 12th node
      //    node12.index = i = 0; parent_of[0] = 11; so we need app[11-1]
      &dat[(parent_of[i]-1) * (*n_char)], // ancestor
      //  To find the number of node 12's children we look in children_of[node12.index]
      //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
      //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
      &dat[(children_of[i + 1 + *n_node] -1) * (*n_char)], // left child
      &dat[(children_of[i + 1] -1) * (*n_char)], // right child
      inapp, end_char
    );
  }
}

void inapp_first_downnode
(int *this, int *left, int *right, 
 int *this_acts, int *l_acts, int *r_acts,
 int *inapp, int *end_char) {
  int i, temp;
  for (i=0; i<*end_char; i++) {
    if ((temp = (left[i] & right[i]))) {
      this[i] = temp;
      if (temp == *inapp) {
        if ((left[i] & ~(*inapp)) && (right[i] & ~(*inapp))) {
          this[i] = (left[i] | right[i]);
        }
      }
    }
    else {
      this[i] = (left[i] | right[i]);
      if ((left[i] & ~(*inapp)) && (right[i] & ~(*inapp))) {
        this[i] = this[i] & ~(*inapp);
      }
    }
    this_acts[i] = (l_acts[i] | r_acts[i]) & ~(*inapp);
//    Rprintf("   - Setting actives (%i+%i) to %i\n", l_acts[i], r_acts[i], this_acts[i]);
  }
}


void inapp_first_downpass
(int *dat, int *act, int *parent, int *child, 
 int *end_char, int *n_char, int *inapp, int *n_edge) {
  int i;  
  for (i = 0; i < *n_edge; i+=2) {
//    Rprintf(" - First downpass at node %i\n", parent[i] - 1);
    inapp_first_downnode(&dat[(parent[i]  - 1) * (*n_char)],
                         &dat[( child[i+1]- 1) * (*n_char)],
                         &dat[( child[i]  - 1) * (*n_char)],
                         &act[(parent[i]  - 1) * (*n_char)],
                         &act[( child[i+1]- 1) * (*n_char)], 
                         &act[( child[i]  - 1) * (*n_char)], inapp, end_char);
  }
}

void inapp_second_downnode
(int *this, int *this_app, int *left, int *right,
 int *this_acts, int *l_acts, int *r_acts,
 int *end_char, int *inapp, int *pars, double *weight, double *w) {
  int i, temp;
  for (i = 0; i < *end_char; ++i) {
    Rprintf("   - this_app = %i, left=%i, right=%i, l_act=%i, r_act=%i\n",
      this_app[i], left[i], right[i], l_acts[i], r_acts[i]);
    if (this_app[i] != *inapp) {
      if ((temp = (left[i] & right[i]))) {
        if (temp & ~(*inapp)) {
          this[i] = temp & ~(*inapp);
        } else {
          this[i] = temp;
        }
      }
      else {
        this[i] = (left[i] | right[i]) & ~(*inapp);
        
        if ((left[i] & ~(*inapp) && right[i] & ~(*inapp))
        ||  (l_acts[i] && r_acts[i])) {
          Rprintf("Addscore %f - %i\n", weight[i], 173);
          (pars[i])++; (*w) += weight[i]; // Add one to tree length
        }
      }
    }
    this_acts[i] = (l_acts[i] | r_acts[i]) & ~(*inapp);
  }
}

void inapp_second_root
(int *dat, int *app, int *inapp, int *end_char) {
  int i;
  for (i=0; i<*end_char; i++) {
    if (app[i] != *inapp) { // Assume applicable if ambiguous.
      dat[i] = dat[i] & ~(*inapp); 
      app[i] = dat[i]; 
    }
  }
}

void inapp_second_downpass
(int *dat, int *app, int *act, int *parent, int *child,
 int *n_edge, int *n_char, int *end_char, int *inapp, 
 int *pars, double *pvec, double *weight) {
  int i, parent_i = 0;
  for (i=0; i<*n_edge; i+=2) {
    parent_i = parent[i];
    pvec[parent_i -1] += pvec[child[i+1]-1] + pvec[child[i]-1];
    Rprintf(" - Calling second downnode at %i:\n", parent_i - 1);
    inapp_second_downnode(
      &dat[(parent_i  -1) * (*n_char)],
      &app[(parent_i  -1) * (*n_char)],
      &dat[(child[i+1]-1) * (*n_char)], 
      &dat[(child[i]  -1) * (*n_char)], 
      &act[(parent_i  -1) * (*n_char)],
      &act[(child[i+1]-1) * (*n_char)], 
      &act[(child[i]  -1) * (*n_char)], 
      end_char, inapp, pars, weight, &pvec[(parent_i-1)]
    );
  }
  inapp_second_root(
    &dat[(parent_i -1) * (*n_char)], 
    &app[(parent_i -1) * (*n_char)], inapp, end_char);
}

void inapp_second_upnode
(int *this, int *ancestor, int *left, int *right, 
 int *this_acts, int *l_acts, int *r_acts,
 int *end_char, int *inapp, int *pars, double *weight, double *w) {
  int i;
  for (i = 0; i < *end_char; ++i) {    
    Rprintf("   - this=%i, anc=%i, left=%i, right=%i, l_act=%i, r_act=%i\n",
            this[i], ancestor[i], left[i], right[i], l_acts[i], r_acts[i]);
    if (this[i] & ~(*inapp)) {
      if (ancestor[i] & ~(*inapp)) {
        if ((ancestor[i] & this[i]) == ancestor[i]) {
          this[i] = ancestor[i] & this[i];
        } else {
          if (left[i] & right[i]) {
            this[i] = (this[i] | (ancestor[i] & left[i] & right[i]));
          }
          else {
            if ((left[i] | right[i]) & (*inapp)) {
              if ((left[i] | right[i]) & ancestor[i]) {
                this[i] = ancestor[i];
              } else {
                this[i] = (left[i] | right[i] | ancestor[i]) & (*inapp);
              }
            } else {
              this[i] = this[i] | ancestor[i];
              if ((ancestor[i] & this[i]) == ancestor[i]) {
                this[i] = ancestor[i] & this[i];
              }
            }
          }
        }
      }
    }
    else {
      if (l_acts[i] && r_acts[i]) {
         Rprintf("Addscore %f - %i\n", weight[i], 251);
        (pars[i])++; (*w) += weight[i]; // Add one to tree length
      }
    }
  }
}

void inapp_second_uppass
(int *dat, int *act, int *parent_of, int *child_of, 
int *end_char, int *n_char, int *n_node, int *inapp, 
int *pars, double *pvec, double *weight) {
  for (int i = 0; i < (*n_node)-1; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // app and app are stored as [0 * n_char] = apps[,1], [1 * n_char] = apps[,2], ....
    Rprintf(" - Calling second upnode at %i, anc=%i\n", parent_of[0] + i, parent_of[i] - 1);
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    inapp_second_upnode(
    //  The position of node 12 in the app array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &dat[(parent_of[0] + i + 1 -1) * (*n_char)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need app[11-1]
    &dat[(parent_of[i]-1) * (*n_char)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &dat[(child_of[i + 1 + *n_node] -1) * (*n_char)], // left child
    &dat[(child_of[i + 1] -1) * (*n_char)], // right child
    
    &act[(parent_of[0] + i + 1 -1) * (*n_char)], // this-actives
    &act[(child_of[i + 1] -1) * (*n_char)], // left child-actives
    &act[(child_of[i + 1 + *n_node] -1) * (*n_char)], // right child-actives
    
    end_char, inapp, pars, weight, &pvec[parent_of[0] + i + 1 -1]);
  }
}

SEXP MORPHYFITCH
(SEXP dat, SEXP nchar, SEXP ntip, 
 SEXP parent, SEXP child, SEXP parent_of, SEXP children_of,
 SEXP weight, SEXP inapp, SEXP inapp_chars) { 
 // Memo: the first 'inapp_chars' characters require the Morphy Treatment.
 // Memo: R plots the first-mentioned child ["left"] below the second-mentioned ["right"]
  double *pvtmp;
  int *data, *appl, *actives, *n_char=INTEGER(nchar), *inappl=INTEGER(inapp), 
      i, n_tips=INTEGER(ntip)[0], *first_applicable=INTEGER(inapp_chars);
  int n_internal = n_tips - 1L, n_edge = (n_tips * 2) - 2L, max_node = n_edge + 1L;
  SEXP RESULT, pars, pscore, DAT, pvec, APPL, ACTIVE;
  PROTECT(RESULT = allocVector(VECSXP, 4));
  PROTECT(pars = allocVector(INTSXP, *n_char));
  PROTECT(pscore = allocVector(REALSXP, 1));
  PROTECT(DAT = allocMatrix(INTSXP, *n_char, max_node));
  PROTECT(pvec = allocVector(REALSXP, max_node));
  PROTECT(APPL = allocMatrix(INTSXP, *n_char, max_node));
  PROTECT(ACTIVE = allocMatrix(INTSXP, *n_char, max_node)); // #TODO One day we should use first_applicable[0] instead
  pvtmp = REAL(pvec);
  for(i=0; i<max_node; i++) pvtmp[i] = 0.0;
  for(i=0; i<*n_char; i++) INTEGER(pars)[i] = 0;
  REAL(pscore)[0] = 0.0;
  data = INTEGER(DAT);
  appl = INTEGER(APPL);
  actives = INTEGER(ACTIVE);
  for(i=0; i<(*n_char * max_node); i++) actives[i] = 0;
  for(i=0; i<(*n_char * n_tips); i++) {
    data[i] = INTEGER(dat)[i];
    // TODO for appl, make i < (*first_applicable * n_tips)
    // The current code is inefficient - we don't need to copy appl for 
    // taxa that are applicable - but then we can copy data more easily. 
    appl[i] = INTEGER(dat)[i];
    actives[i] = INTEGER(dat)[i];
  }

  inapp_first_downpass(data, actives, INTEGER(parent), INTEGER(child), 
                       first_applicable, n_char, inappl, &n_edge);
  inapp_first_uppass(data, appl, INTEGER(parent_of), INTEGER(children_of),
                     first_applicable, n_char, &n_internal, inappl);
  // TODO!: inapp_update_tips   (appl, first_applicable, INTEGER(parent), INTEGER(child), INTEGER(n_edge)); 
  inapp_second_downpass(data, appl, actives, INTEGER(parent), INTEGER(child), 
                        &n_edge, n_char, first_applicable, inappl, 
                        INTEGER(pars), pvtmp, REAL(weight));
  inapp_second_uppass  (data, actives, INTEGER(parent_of), INTEGER(children_of), 
                        first_applicable, n_char, &n_internal, inappl, 
                        INTEGER(pars), pvtmp, REAL(weight));

  app_fitch_downpass(data, INTEGER(parent), INTEGER(child),
                     first_applicable, n_char, &n_edge, 
                     inappl, INTEGER(pars), REAL(weight), pvtmp); // No need for an up-pass: all scoring on way down.
  
  *(REAL(pscore)) = pvtmp[n_tips];
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(7);
  return(RESULT);
}