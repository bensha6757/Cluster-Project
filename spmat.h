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
	/* An array s.t. spmatSize[i] == the count of non-zero elements in row i */ 
	int_vector spmatSize;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const double *row, int i);

	/*Gets row i from the matrix */
	/*void (*get_row)(const struct _spmat *A, int i, int_vector K, vector row, num M);*/
	void (*get_row)(const struct _spmat *A, int i, vector row);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	struct _spmat* (*create_sub_mat)(const struct _spmat  *A, Subgroup g, int sizeG);
	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *private;
	
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n);

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

/** Set e_i to be a n-length vector s.t. 
 * e[j]==1 if i==j, else e[j]==0. */
void get_basis_unit_vec(vector *e_i, num i, num n);

num sum_array(int_vector vec, num len);

#endif /* SPMAT_H_ */
