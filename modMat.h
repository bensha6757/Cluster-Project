
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

#define USE_SPMAT_LINKED TRUE
#define MODULARITY_INIT -1

/*
 * A struct representing Modularity Matrix B[g].
 */
typedef struct _modmat {
	spmat *A; 				/* Network adjacency matrix in form of spmat */
	int_vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	num M;					/* Total sum of degrees in the network*/
	num currM;				/* Total sum of degrees in sub-network defined by some Subgroup g */
	num gSize;	 			/* Size of matrix, for any Subgroup g it is reduced to */
	double one_norm;		/* The 1-norm of the matrix, i.e. max_i(sum_j(|B_ij|)) */ 

	/*free all resources used by a ModMat instance */
	void (*free)(struct _modmat *B);

	/* Multiply ModMat with vector v and store result.
	 */
	void (*mult)(const struct _modmat *B, const double *v, double *result, boolean shift);

	/** Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
	*	the functions can be called from the Maximization algorithm, where only one entry of s is being changed, so the spmat module can calculate it smartly and quickly.
	* 	@param move_vertex - if MODULARITY_INIT, compute modularity w.r.t to s. Otherwise, move vertex s[move_vertex] temporarily and compute.
	*/
	double (*get_modularity)(struct _modmat *B, vector s, int move_vertex);

} modMat;

/** Compute dot-product of two vectors of real numbers of size d.
 */
double dot_prod(vector v, vector u, num d);

/* Allocate a new, empty instance of Modularity Matrix */
modMat* allocate_mod_mat(num n, num nnz, boolean isSub);

/* Constructor, creates a new sub matrix B_hat[g], based on B and g */
modMat *create_Sub_Matrix(modMat *B, Subgroup g, num sizeG);

/* Computes and sets 1-norm of B, i.e. max_i(sum_j(abs(B_ij))) */
void set_1_norm(modMat *B);

#endif /* MODMAT_H_ */
