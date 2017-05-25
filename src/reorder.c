#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void order_edges(int *parent, int *child, const int *n_node, const int *n_edge) {
  int i, q_pos = 0, o_node, root_node;
  int * start_p = calloc(*n_edge, sizeof(int));
  int * start_c = calloc(*n_edge, sizeof(int));
  int * child_l  = calloc(*n_node, sizeof(int));
  int * child_r  = calloc(*n_node, sizeof(int));
  int * queue_p  = calloc(*n_node, sizeof(int));
  int * queue_c  = calloc(*n_node, sizeof(int));
  // TODO check that calloc has returned a non-null pointer; clean-up and exit if calloc has failed
  root_node = *n_node + 2L;
  for (i = 0; i < *n_edge; i++) {
    // Initialize
    start_p[i] = parent[i];
    start_c[i] = child[i];
    q_pos = parent[i] - root_node;
    if (child_l[q_pos]) {
      if (child_l[q_pos] < child[i] && child[i] > root_node) {
        child_r[q_pos] = child[i];
      } else {
        child_r[q_pos] = child_l[q_pos];
        child_l[q_pos] = child[i];
      }
    } else {
      child_l[q_pos] = child[i];
    }
  }
  o_node = root_node;
  q_pos = 0;
  for (i = 0; i < *n_edge; i++) {
    if (o_node < root_node) { // We've just reached a tip
      parent[i] = queue_p[--q_pos];
      child[i] = queue_c[q_pos];
      o_node = child[i];
    } else { // We're at an internal node
      parent[i] = o_node;
      child[i]  = child_l[o_node - root_node];
      queue_p[q_pos] = o_node;
      queue_c[q_pos++] = child_r[o_node - root_node];
      o_node = child_l[o_node - root_node];
    }
  }
  free(start_p);
  free(start_c);
  free(child_l);
  free(child_r);
  free(queue_p);
  free(queue_c);
}

void number_nodes(int *parent, int *child, const int *root_node, const int *n_edge) {
  int i, next_node, n_allnodes = *n_edge + 1L;
  int * renumber = calloc(n_allnodes, sizeof(int));
  // TODO check that calloc has returned a non-null pointer; clean-up and exit if calloc has failed
  next_node = *root_node;
  for (i = 0; i < n_allnodes; i++) renumber[i] = i + 1L;
  for (i = 0; i < *n_edge; i++) {
    if (child[i] > *root_node) renumber[child[i]-1L] = ++(next_node);
  }
  for (i = 0; i < *n_edge; i++) {
    parent[i] = renumber[parent[i]-1L];
    child[i] = renumber[child[i]-1L];
  }
  free(renumber);
}

void tabulate
(const int *x, const int *n, const int *nbin, int *ans) {
    int i, tmp;
    for (i=0; i < *nbin; i++) ans[i]=0L; 
    for (i=0; i < *n; i++) {
        tmp = x[i];
        if( (tmp>0) & (tmp<(*nbin+1L)) )   
        ans[tmp-1L] ++;
    }
}

void phangorn_reorder
(const int *from, const int *to, const int *n, 
 const int *sumNode, int *neworder, int *root) { 
    int i, j, sum=0, k, Nnode, ind, *ord, *csum, *tips, *stack, z=0;  // l, 
    double *parent;
    int m=sumNode[0];
    parent = (double *) R_alloc((*n), sizeof(double));
    tips = (int *) R_alloc(m, sizeof(int));
    ord = (int *) R_alloc((*n), sizeof(int));
    csum = (int *) R_alloc( (m+1), sizeof(int));
    stack = (int *) R_alloc(m, sizeof(int));
    for(j=0;j<(*n);j++) parent[j] = (double)from[j];
   
    for(j=0;j<(*n);j++) ord[j] = j;
    for(j=0;j<m;j++) tips[j] = 0;
        
    rsort_with_index(parent, ord, *n);
    tabulate(from, n, sumNode, tips);
    csum[0]=0;
    for(i=0;i<(*sumNode);i++){
        sum+=tips[i];                 
        csum[i+1] = sum;                        
    }      
    k = (*n)-1;
    Nnode = 0;
    stack[0] = *root;
    
    while(z > -1){
        j=stack[z];          
        if(tips[j]>0){   
            for(i=csum[j];i<csum[j+1];i++){
                ind = ord[i];                     
                neworder[k] = ind + 1;        
                stack[z] = to[ind]-1;
                k -=1;
                z++;                            
            }         
            Nnode += 1; 
            }
        z--;       
    }                
    root[0]=Nnode;     
}