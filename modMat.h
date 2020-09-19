
#ifndef MODMAT_H_
#define MODMAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "Spmat.h"
#include "IO_Mem_Errors.h"
#include "Types.h"

/*A choice of implementation for spmat module. Can be changed */
#define USE_SPMAT_LINKED (TRUE)

#define MODULARITY_INIT (-1)
#define SHIFT (TRUE)
#define NO_SHIFT (FALSE)

/** A struct representing Modularity Matrix B. */
typedef struct _modmat {
	spmat *A; 				/* Network adjacency matrix in form of spmat */
	int_vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	num M;					/* Total sum of degrees in the network*/
	num currM;				/* Total sum of degrees in sub-network defined by some Subgroup g */
	num gSize;	 			/* Size of matrix, for any Subgroup g it is reduced to */
	double one_norm;		/* The 1-norm of the matrix, i.e. max_i(sum_j(|B_ij|)) */ 

	/** Free all resources used by a ModMat instance. */
	void (*free)(struct _modmat *B);

	/** Multiply ModMat with vector v and store output in result.
	 *  @param shift - a flag that can be assigned SHIFT or NO_SHIFT. If SHIFT - The modMat instance multiplied will be shifted by one_norm value.
	 */
	void (*mult)(const struct _modmat *B, const double *v, double *result, boolean shift);

	/** Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
	*	the functions can be called from the Maximization algorithm, where only one entry of s is being changed, so the spmat module can calculate it smartly and quickly.
	* 	@param move_vertex - if MODULARITY_INIT, compute modularity w.r.t to s. Otherwise, move vertex s[move_vertex] temporarily and compute.
	*/
	double (*get_modularity)(struct _modmat *B, vector s, int move_vertex);

} modMat;

/** Compute dot-product of two vectors of real numbers of dimension d.
 */
double dot_prod(vector v, vector u, num d);

/** Allocate a new, empty instance of Modularity Matrix.
 *  @param n - the dimension of the modularity matrix to be allocated (size n by n).
 * 	@param nnz - the number of non-zero elements in the mod. matrix. Optional if isSub == TRUE.
 *  @param isSub - a boolean. If TRUE, a the field spmat A is assumed to be allocated and assigned externally.
 */
modMat* allocate_mod_mat(num n, num nnz, boolean isSub);

/** Constructor, creates a new sub matrix B[g], based on B and g, 
 *  a subgroup of indices in increasing order, of size sizeG. */
modMat *create_Sub_Matrix(modMat *B, Subgroup g, num sizeG);

/** Computes and sets 1-norm of B, i.e. max_i(sum_j(abs(B_ij))) */
void set_1_norm(modMat *B);

#endif /* MODMAT_H_ */
