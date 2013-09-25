#define USE_RINTERNALS
#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

void order_edges(int *parent, int *child, int *n_node, int *n_edge) {
  int i, q_pos = 0, o_node, root_node;
  int * start_p = calloc(*n_edge, sizeof(int));
  int * start_c = calloc(*n_edge, sizeof(int));
  int * child_l  = calloc(*n_node, sizeof(int));
  int * child_r  = calloc(*n_node, sizeof(int));
  int * queue_p  = calloc(*n_node, sizeof(int));
  int * queue_c  = calloc(*n_node, sizeof(int));
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
}

void number_nodes(int *parent, int *child, int *root_node, int *n_edge) {
  int i, next_node, n_allnodes = *n_edge + 1L;
  int * renumber = calloc(n_allnodes, sizeof(int));
  next_node = *root_node;
  for (i = 0; i < n_allnodes; i++) renumber[i] = i + 1L;
  for (i = 0; i < *n_edge; i++) {
    if (child[i] > *root_node) renumber[child[i]-1L] = ++(next_node);
  }
  for (i = 0; i < *n_edge; i++) {
    parent[i] = renumber[parent[i]-1L];
    child[i] = renumber[child[i]-1L];
  }
}