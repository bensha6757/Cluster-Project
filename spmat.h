#ifndef SPMAT_H_
#define SPMAT_H_

#include "io_mem_errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "Types.h"

typedef struct _spmat {
	/* Matrix size (n*n) */
	int n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const double *row, int i);

	/*Gets row i from the matrix */
	void (*get_row)(const struct _spmat *A, int i, double *row);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *private;
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

/** Returns a populated sparse matrix (linked-list impl.) of size n 
 * with rows reduced to elements according to col. indices in Subgroup g. */
spmat* create_sub_sparse_matrix_linked(spmat *A, Subgroup g, int n , int_vector spmatSize);

/** Returns a populated sparse matrix (arrays impl.) of size n 
 * with rows reduced to elements according to col. indices in Subgroup g. */
spmat* create_sub_sparse_matrix_array(spmat *A, Subgroup g, int sizeG, int_vector spmatSize);

/** Set e_i to be a n-length vector s.t. 
 * e[j]==1 if i==j, else e[j]==0. */
void get_basis_unit_vec(vector *e_i, num i, num n);

#endif /* SPMAT_H_ */
