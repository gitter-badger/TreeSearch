void report_calloc_error() {
  Rprintf("Error allocating memory in calloc");  
}

static R_NativePrimitiveArgType order_edges_number_nodes_t[] = {
  INTSXP, INTSXP, INTSXP
};
extern void order_edges_number_nodes(int *parent, int *child, const int *n_edge)
{
  int i, queue_pos = 0, o_node, next_node;
  const int n_node = *n_edge / 2;
  const int n_allnodes = *n_edge + 1L, root_node = n_node + 2L;
  int * start_p = calloc(*n_edge, sizeof(int)), /* calloc zero-initializes */
      * start_c = calloc(*n_edge, sizeof(int)),
      * child_l = calloc( n_node, sizeof(int)),
      * child_r = calloc( n_node, sizeof(int)),
      * queue_p = calloc( n_node, sizeof(int)),
      * queue_c = calloc( n_node, sizeof(int)),
      * renumber = calloc(n_allnodes, sizeof(int));
  if (queue_c == NULL) report_calloc_error(); 
  /* If queue_c has calloc'ed, presume other callocs were successful too */
    
  for (i = 0; i < *n_edge; i++) {
    /* Initialize */
    start_p[i] = parent[i];
    start_c[i] = child[i];
    queue_pos = parent[i] - root_node;
    if (child_l[queue_pos]) {
      if (child_l[queue_pos] < child[i] && child[i] > root_node) {
        child_r[queue_pos] = child[i];
      } else {
        child_r[queue_pos] = child_l[queue_pos];
        child_l[queue_pos] = child[i];
      }
    } else {
      child_l[queue_pos] = child[i];
    }
  }
  o_node = root_node;
  queue_pos = 0;
  for (i = 0; i < *n_edge; i++) {
    if (o_node < root_node) { /* We've just reached a tip */
      parent[i] = queue_p[--queue_pos];
      child[i] = queue_c[queue_pos];
      o_node = child[i];
    } else { /* We're at an internal node */
      parent[i] = o_node;
      child[i]  = child_l[o_node - root_node];
      queue_p[queue_pos] = o_node;
      queue_c[queue_pos++] = child_r[o_node - root_node];
      o_node = child_l[o_node - root_node];
    }
  }
  free(start_p);
  free(start_c);
  free(child_l);
  free(child_r);
  free(queue_p);
  free(queue_c);
  
  /* Now number nodes: */
  if (renumber != NULL) {
    next_node = root_node;
    for (i = 0; i < n_allnodes; i++) renumber[i] = i + 1;
    for (i = 0; i < *n_edge; i++) {
      if (child[i] > root_node) renumber[child[i]-1] = ++(next_node);
    }
    for (i = 0; i < *n_edge; i++) {
      parent[i] = renumber[parent[i]-1L];
      child[i] = renumber[child[i]-1L];
    }
    free(renumber);
  } else {
    report_calloc_error();
  }
}

extern SEXP RENUMBER_TREE(SEXP parent, SEXP child, SEXP ned) {
  int i;
  const int n_edge = INTEGER(ned)[0];
  SEXP RESULT;
  PROTECT(RESULT = allocVector(INTSXP, n_edge * 2));
  for (i = 0; i < n_edge; i++) {
    INTEGER(RESULT)[i] = INTEGER(parent)[i];
    INTEGER(RESULT)[i + n_edge] = INTEGER(child)[i];
  }
  
  order_edges_number_nodes(INTEGER(RESULT), INTEGER(RESULT) + n_edge, &n_edge);
  
  UNPROTECT(1);
  return(RESULT);
}

extern SEXP RENUMBER_EDGES(SEXP parent, SEXP child, SEXP ned) {
  int i;
  const int n_edge = INTEGER(ned)[0];
  SEXP RESULT, PARENT, CHILD;
  PROTECT(RESULT = allocVector(VECSXP, 2L));
  PROTECT(PARENT = allocVector(INTSXP, n_edge));
  PROTECT(CHILD  = allocVector(INTSXP, n_edge));
  for (i = 0; i < n_edge; i++) {
    INTEGER(PARENT)[i] = INTEGER(parent)[i];
    INTEGER(CHILD )[i] = INTEGER(child )[i];
  }
  
  order_edges_number_nodes(INTEGER(PARENT), INTEGER(CHILD), &n_edge);
    
  SET_VECTOR_ELT(RESULT, 0, PARENT);
  SET_VECTOR_ELT(RESULT, 1, CHILD);
  UNPROTECT(3);
  return(RESULT);
}
