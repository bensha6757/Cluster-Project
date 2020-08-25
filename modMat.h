/*
 * modMat.h
 *
 *  Created on: 22 ����� 2020
 *      Author: ��
 */

#ifndef MODMAT_H_
#define MODMAT_H_

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

typedef size_t* int_vector;
typedef int_vector Subgroup;

/*
 * A struct representing Modularity Matrix B[g].
 */
typedef struct _modmat {
	spmat *A; 			/* Network adjacency matrix in form of spmat */
	int_vector spmatSize; 	/* a vector of the spmat rows sizes for internal use */
	int_vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	size_t M; 			/* Total sum of degrees in the network*/
	size_t gSize;	 	/* Size of subgroup g  */
	Subgroup g; 		/* Array of relevant indices of submatrix */

	/*free all resources used by ModMat instance */
	void (*free)(struct _modmat *B, Subgroup g);

	/*Get row i of B*/
	void (*get_row)(const struct _modmat *B, size_t i, double *row);

	void (*mult)(const struct _modmat *B, const struct _modmat *Bg, const double *v, double *result);
	/*multiply ModMat with vector v and store result.
	 * Can be extended for submatrix B[g].
	 */

} modMat;

/* Allocate a new, empty instance of Modularity Matrix */
modMat *allocate_mod_mat(int n);

/* Get a symmetric B_hat instance's 1-norm, i.e. ||B_hat||_1 = max_j(sum_i(|B_hat_ij|)).
 */
double get_1_norm(modMat *B);

/* Read the number n=|V| of network from open file to memory*/
size_t read_totalV_from_file(FILE *input);

/* Load an open input file into previously allocated ModMat B*/
void load_mod_matrix_from_file(FILE *input, modMat *B);

#endif /* MODMAT_H_ */
