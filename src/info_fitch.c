#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>
/*

double factorialScaledBy(double m, unsigned int n) {
  if (n < 2) return m*1;
  return factorialScaledBy(n*m, n-2);
}
double dfactorial (unsigned int n) {
  return factorialScaledBy(1, n);
}

double dfact (unsigned int n) {
  if (n == 0 || n == 1) return 1;
  double fact = 1.0;
  for (unsigned int i = n; i > 1; i -= 2) fact = fact * i;
  return fact;
}
double n_unrooted (unsigned int n) {return dfact(2 * n - 5);} // n > 1
double n_rooted   (unsigned int n) {return dfact(2 * n - 3);} // n > 1

double info2(int split1, int split2) {
  double poss_trees = n_unrooted(split1 + split2);
  double valid_trees = n_rooted(split1) + n_rooted(split2);
  return -log2(valid_trees/poss_trees);
}

void info_content(int *zone, int *n_rows, int *n_tip, int *n_node, int *max_node, double *inf) {
  int n_splits = 0;
  int arr[*n_node];
  for (int i = 0; i < *n_node; i++) arr[i] = 0;
  for (int i = 0; i < *n_rows; i++) {
    n_splits = 0;
    for (int j = *n_tip; j < *max_node; j++) {
      if (zone[j * i]) {
        arr[n_splits++] = zone[j * i];
      }
    }
    switch (n_splits) {
      case 0: break;
      case 1: break;
      case 2: inf[i] = info2(arr[0], arr[1]); break;
      default: inf[i] = 999;
    }
  }
}
*/

void fitch_inf_upnode(int *down, int *up, int *ancestor, int *child_q, int *child_r, int *n_rows) {
  for (int k = 0; k < (*n_rows); k++) {
    if ((down[k] & ancestor[k]) == ancestor[k]) {
      up[k] = down[k] & ancestor[k];
    } else if (child_q[k] & child_r[k]) {
      up[k] = down[k] & ancestor[k];
    } else {
      up[k] = down[k] | (ancestor[k] & (child_q[k] | child_r[k]));
    }
  }
}

void fitch_inf_uproot(int *down, int *up, int *n_rows) {
  for (int k = 0; k < (*n_rows); k++) up[k] = down[k];
}

void fitch_inf_uppass(int *dat, int *up, int *app, int *n_rows, int *parent_of, int *children_of, int *n_node, int *inapp) {
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  fitch_inf_uproot(&dat[(parent_of[0]-1) * (*n_rows)], &up[(parent_of[0]-1) * (*n_rows)], n_rows);
  for (int i = 0; i < (*n_node)-1; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // app and app are stored as [0 * n_rows] = apps[,1], [1 * n_rows] = apps[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    fitch_inf_upnode(
    //  The position of node 12 in the app array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &dat[(parent_of[0] + i + 1 -1) * (*n_rows)], // this_start
    &up [(parent_of[0] + i + 1 -1) * (*n_rows)], // this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need app[11-1]
    &dat[(parent_of[i]-1) * (*n_rows)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &dat[(children_of[i + 1]-1) * (*n_rows)], // child q
    &dat[(children_of[i + 1 + *n_node]-1) * (*n_rows)], // child r
    n_rows);
  }
}

void fitch_inf_downnode(int *app, int *this, int *child_q, int *child_r, int *n_rows, int *inapp) {
  int tmp;
  for (int k = 0; k < (*n_rows); k++) {
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
    }
  }
}

void fitch_inf_downpass(int *dat, int *app, int *n_rows, int *parent, int *child, int *n_edge, int *inapp) {
  int parent_i = parent[0];
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    fitch_inf_downnode(&app[(parent_i -1)*(*n_rows)], &dat[(parent_i -1)* (*n_rows)], &dat[(child[i]-1) * (*n_rows)], &dat[(child[i+1]-1) * (*n_rows)], n_rows, inapp);
  }
}

void inf_upnode(int *this, int *ancestor, int *child_q, int *child_r, int *n_rows) {
  for (int k = 0; k < (*n_rows); k++) this[k] = ((child_q[k] + child_r[k] + ancestor[k] >= 0) ? 1 : -1);
}

void inf_uproot(int *this, int *n_rows) {
  for (int k = 0; k < (*n_rows); k++) this[k] = ((this[k] > -1) ? 1 : -1); // Assume applicable if ambiguous.
}

void inf_uppass(int *app, int *n_rows, int *parent_of, int *children_of, int *n_node) {
  inf_uproot(&app[(parent_of[0]-1) * (*n_rows)], n_rows);
  // parent_of's first member is a child of the root node.  Thus parent_of[0] = root.number
  for (int i = 0; i < (*n_node)-1; i++) {
    // parent_of is stored as 1L, 1R, 2L, 2R, 3L, 3R, 4L, 4R, ... nL, nR.  (The root node has no parent.)
    // children_of is stored as 0L, 1L, 2L, ... nL, 0R, 1R, 2R, 3R, ..., nR
    // app and app are stored as [0 * n_rows] = apps[,1], [1 * n_rows] = apps[,2], ....
    
    // Worked example assumes that root node = 11 and i = 0, meaning 'look at node 12' [the first of 11's children].
    inf_upnode(
    //  The position of node 12 in the app array is:
    //    root.number (counting from 1) + i = i.e. position[11], the 12th position
    &app[(parent_of[0] + i + 1 -1) * (*n_rows)], // this_start, will become this_finish
    //  To find the number of node 12's parent we look in parent_of[node12.index]
    //    parent_of[0] is the parent of node [root + i] = 12th node
    //    node12.index = i = 0; parent_of[0] = 11; so we need app[11-1]
    &app[(parent_of[i]-1) * (*n_rows)], // ancestor
    //  To find the number of node 12's children we look in children_of[node12.index]
    //    children_of[0, *n_node + 0] are the two children of node [root + i] = 12
    //    node12.index = i = 0; children_of[0*2] = Q; children_of[0*2 + 1] = R
    &app[(children_of[i + 1]-1) * (*n_rows)], // child q
    &app[(children_of[i + 1 + *n_node]-1) * (*n_rows)], // child r
    n_rows);
  }
}

void inf_downnode(int *dat, int *child_q, int *child_r, int *n_rows) {
  int tmp;
  for (int k = 0; k < (*n_rows); k++) {
    dat[k] = ((tmp = (child_q[k] + child_r[k])) ? (tmp > 0 ? 1 : -1) : 0);
  }
}

void inf_downpass(int *app, int *n_rows, int *parent, int *child, int *n_edge) {
  for (int i = 0; i < *n_edge; i+=2) inf_downnode(&app[(parent[i]-1) * (*n_rows)], &app[(child[i]-1) * (*n_rows)], &app[(child[i+1]-1) * (*n_rows)], n_rows);
}

void change_edge(int *parent_up, int *child_down, int *child_up, int *n_rows, int *changes) {
  // Algorithm 3 in http://phylo.bio.ku.edu/slides/BIOL848-lec7-ML-MPModelAlgorithms.pdf
  unsigned int s, t;
  for (int k = 0; k < (*n_rows); k++) {
    for (s = 1<<31; s != 0; s >>= 1) {
      if (parent_up[k] & s) { // for each state s in the up-pass set of p
        if (child_down[k] & s) { //if s is in the down-pass set of c
          // s->s is a state to state transition across branch e
        } else {
          for (t = 1<<31; t != 0; t >>= 1) { // for each state t in the down-pass set of c
            changes[k] |= (child_down[k] & t);
          }
          if (child_up[k] & s) {
            // s->s is a state to state transition across branch e
          }
        }
      }
    }
  }
}

void changes (int *downpass, int *uppass, int *n_rows, int *parent, int *child, int *n_edge, int *changes) {
  for (int i = 0; i < *n_edge; i++) change_edge(
    &uppass[(parent[i]-1) * (*n_rows)],
    &downpass[(child[i]-1) * (*n_rows)],
    &uppass[(child[i]-1) * (*n_rows)],
    n_rows, &changes[i * (*n_rows)]);
}

void downzone (int *this, int *child_q, int *child_r, int *change_q, int *change_r, int *up_q, int *up_r, int *n_rows, int *inapp) {
  for (int k = 0; k < (*n_rows); k++) {
    this[k] = 0;
    if (!change_q[k]) {
      if (up_q[k] != *inapp) this[k] += child_q[k];
      child_q[k] = 0;
    }
    if (!change_r[k]) {
      if (up_r[k] != *inapp) this[k] += child_r[k];
      child_r[k] = 0;
    }
  }
}

void zones (int *uppass, int *change, int *n_rows, int *parent, int *child, int *child_edge, int *n_edge, int *n_node, int *inapp, int *zone) {
  int parent_i = parent[0];
  int child_q = child[0];
  int child_r = child[0+1];
  for (int i = 0; i < *n_edge; i+=2) {
    parent_i = parent[i];
    child_q = child[i];
    child_r = child[i+1];
    downzone(
      &zone[(parent_i -1)* (*n_rows)], &zone[(child_q-1) * (*n_rows)], &zone[(child_r-1) * (*n_rows)],
      &change[i * (*n_rows)], &change[(i+1) * (*n_rows)], 
      &uppass[(child_q-1) * (*n_rows)], &uppass[(child_r-1) * (*n_rows)],
      n_rows, inapp
    );
  }
}


SEXP FITCHTRANS(SEXP dat, SEXP nrx, SEXP parent, SEXP child, SEXP parent_of, SEXP children_of, SEXP child_edge, SEXP n_edge, SEXP n_node, SEXP max_node, SEXP n_tip, SEXP inapp) {   
  int *data, *uppass, *appl, *change, *zone, *n_rows=INTEGER(nrx), e=INTEGER(n_edge)[0], *inappl=INTEGER(inapp),  m=INTEGER(max_node)[0], i, n=INTEGER(n_tip)[0];
//  double *inf;
  SEXP RESULT, DAT, UPPASS, APPL, CHANGE, ZONES;
//SEXP INF;
  PROTECT(RESULT = allocVector(VECSXP, 5));
  PROTECT(DAT = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(UPPASS = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(APPL = allocMatrix(INTSXP, n_rows[0], m));
  PROTECT(CHANGE = allocMatrix(INTSXP, n_rows[0], e));
  PROTECT(ZONES = allocMatrix(INTSXP, n_rows[0], m));
//  PROTECT(INF = allocVector(REALSXP, n_rows[0]));
  data = INTEGER(DAT);
  uppass = INTEGER(UPPASS);
  appl = INTEGER(APPL);
  change = INTEGER(CHANGE);
  zone  = INTEGER(ZONES);
//  inf   = REAL(INF);
  for (i=0; i<(*n_rows * n); i++) {
    data[i] = INTEGER(dat)[i];
    uppass[i] = INTEGER(dat)[i];
    appl[i] = ((data[i] & ~*inappl) ? 1 : 0) - ((data[i] & *inappl) ? 1 : 0);
    zone[i] = 1;
  }
  for (i=0; i<(*n_rows * e); i++) {
    change[i] = 0;
  }
//  for (i = 0; i < *n_rows; i++) {
//    inf[i] = 0;
//  }
  inf_downpass(appl, n_rows, INTEGER(parent), INTEGER(child), INTEGER(n_edge));
  inf_uppass(appl, n_rows, INTEGER(parent_of), INTEGER(children_of), INTEGER(n_node));
  fitch_inf_downpass(data, appl, n_rows,  INTEGER(parent), INTEGER(child), INTEGER(n_edge), INTEGER(inapp));
  fitch_inf_uppass(data, uppass, appl, n_rows, INTEGER(parent_of), INTEGER(children_of), INTEGER(n_node), INTEGER(inapp));
  changes(data, uppass, n_rows, INTEGER(parent), INTEGER(child), INTEGER(n_edge), change);
  zones(uppass, change, n_rows, INTEGER(parent), INTEGER(child), INTEGER(child_edge), INTEGER(n_edge), INTEGER(n_node), INTEGER(inapp), zone);
  //info_content(zone, n_rows, INTEGER(n_tip), INTEGER(n_node), INTEGER(max_node), inf);
  
  SET_VECTOR_ELT(RESULT, 0, DAT);
  SET_VECTOR_ELT(RESULT, 1, UPPASS);
  SET_VECTOR_ELT(RESULT, 2, APPL);
  SET_VECTOR_ELT(RESULT, 3, CHANGE);
  SET_VECTOR_ELT(RESULT, 4, ZONES);
 // SET_VECTOR_ELT(RESULT, 5, INF);
  UNPROTECT(6);
  return(RESULT);
}