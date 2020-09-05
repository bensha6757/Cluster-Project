
#ifndef MODMAT_H_
#define MODMAT_H_

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "io_mem_errors.h"
#include <math.h>
#include "Types.h"

#define USE_SPMAT_LINKED 0


/*
 * A struct representing Modularity Matrix B[g].
 */
typedef struct _modmat {
	spmat *A; 				/* Network adjacency matrix in form of spmat */
	int_vector spmatSize; 	/* a vector of the spmat rows sizes for internal use */
	int_vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	num M;					/* Total sum of degrees in the network*/
	num currM;				/* Total sum of degrees in sub-network defined by some Subgroup g */
	num gSize;	 			/* Size of matrix, for any Subgroup g it is reduced to */
	double one_norm;		/* The 1-norm of the matrix, i.e. max_i(sum_j(|B_ij|)) */ 

	/*free all resources used by ModMat instance */
	void (*free)(struct _modmat *B);

	/*Get row i of B*/
	void (*get_row)(const struct _modmat *B, num i, double *row);

	/* Multiply ModMat with vector v and store result.
	 */
	void (*mult)(const struct _modmat *B, const double *v, double *result, boolean shift);


} modMat;

/* Allocate a new, empty instance of Modularity Matrix */
modMat *allocate_mod_mat(num n, num m);

/* constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */
modMat *create_Sub_Matrix(modMat *B, Subgroup g, num sizeG);

/* Computes and sets 1-norm of B, i.e. max_i(sum_j(abs(B_ij))) */
void set_1_norm(modMat *B);

#endif /* MODMAT_H_ */
