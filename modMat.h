/*
 * modMat.h
 *
 *  Created on: 22 באוג׳ 2020
 *      Author: גל
 */

#ifndef MODMAT_H_
#define MODMAT_H_

#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

typedef size_t* vector;
typedef vector Subgroup;

/*
 * A struct representing Modularity Matrix B[g].
 */
typedef struct _modmat {
	spmat *A; 			/* Network adjacency matrix in form of spmat */
	vector spmatSize; 	/* a vector of the spmat rows sizes for internal use */
	vector K; 			/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	size_t M; 			/* Total sum of degrees in the network*/
	size_t gSize;	 	/* Size of subgroup g  */
	Subgroup g; 		/* Array of relevant indices of submatrix */


	void (*free)(struct _modmat *B, Subgroup g); /*free all resources used by ModMat instance */

	void (*mult)(const struct _modmat *B, const struct _modmat *Bg, const double *v, double *result);
	/*multiply ModMat with vector v and store result.
	 * Can be extended for submatrix B[g].
	 */

} modMat;


modMat *allocateModMat(int n);

/*Copy row i of B_hat matrix of size B->gSize to row vector */
void getBhatRow(modMat *B, size_t i, double *row);

/* Get a symmetric B_hat instance's 1-norm, i.e. ||B_hat||_1 = max_j(sum_i(|B_hat_ij|)).
 */
double getOneNorm(modMat *B);

/* Read the number n=|V| of network from open file to memory*/
void readVnumFromFile(FILE *input, size_t *n);

/* Load an open input file into previously allocated ModMat B*/
void loadModMatrixFromFile(FILE *input, modMat *B);

#endif /* MODMAT_H_ */
