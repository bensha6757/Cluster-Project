
#ifndef MODMAT_H_
#define MODMAT_H_

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "io_mem_errors.h"
#include <math.h>
#include "Types.h"

#define USE_LINKED 1


/*
 * A struct representing Modularity Matrix B[g].
 */
typedef struct _modmat {
	spmat *A; 				/* Network adjacency matrix in form of spmat */
	int_vector spmatSize; 	/* a vector of the spmat rows sizes for internal use */
	int_vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	size_t M;				/* Total sum of degrees in the network*/
	size_t gSize;	 		/* Size of subgroup g  */
	Subgroup g; 			/* Array of relevant indices of submatrix */
	double one_norm;		/* The 1-norm of the matrix, i.e. max_i(sum_j(|B_ij|)) */ 

	/*free all resources used by ModMat instance */
	void (*free)(struct _modmat *B);

	/*Get row i of B*/
	void (*get_row)(const struct _modmat *B, size_t i, double *row);

	/* Multiply ModMat with vector v and store result.
	 * Can be extended for submatrix B[g].
	 */
	void (*mult)(const struct _modmat *B, const double *v, double *result);


} modMat;

/* Allocate a new, empty instance of Modularity Matrix */
modMat *allocate_mod_mat(int n);

/* constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */
modMat *create_Sub_Matrix(modMat *B, Subgroup g, size_t sizeG, int impl_flag);

/* Implements multiplication of B_hat[g] with a vector by
 * using several mult. functions and adding results together */
void mult_B_hat_g(modMat *B, modMat *Bg, double *v, double *result);

/* Computes and sets 1-norm of B, i.e. max_i(sum_j(abs(B_ij))) */
void set_1_norm(modMat *B);

#endif /* MODMAT_H_ */
