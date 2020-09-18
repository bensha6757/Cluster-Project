#ifndef SPMAT_H_
#define SPMAT_H_

#include "IO_Mem_Errors.h"
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include "Types.h"

/************************************************************************************************************
 * A data structure representing a sparse matrix, allowing fast creation (adding rows) and multiplication by
 * a vector of real numbers.                              
 * Additional methods leverage internal implementation to enhance computations associated with network 
 * modularity algorithms.
 ************************************************************************************************************/
typedef struct _spmat {
	/** Matrix dimension (n*n) */
	int n;

	/** An array s.t. spmatSize[i] == the number of non-zero elements in row i */
	int_vector spmatSize;

	/** Adds row i the matrix. Called before any other call,
	 * 	exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const double *row, int i);

	/** Gets row i from the matrix */
	void (*get_row)(const struct _spmat *A, int i, vector row);

	/** Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/** Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	/** Creates a submatrix of A of dim == sizeG, which consists only of rows and columns whose index in Subgroup g */
	struct _spmat* (*create_sub_mat)(const struct _spmat  *A, Subgroup g, int sizeG);
	
	/** Computes the score (modularity difference) of a sub-network w.r.t to a 2-division vector s and a temporarily moved vertex.
	 * 	Leverages spmat internal representation for fast computation.
	 * 	Assumes all non-zero elements of spmat are 1's.
	 * 	@param K - The parent network's degrees array.
	 * 	@param M - The parent network's sum of degrees.
	 */ 
	double (*get_modularity_score)(const struct _spmat *A, vector s, int moved_v, int_vector K, num M);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *private;
	
} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* spmat_allocate_list(int n);
/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz);

#endif /* SPMAT_H_ */
