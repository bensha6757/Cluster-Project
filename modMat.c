/*
 * modMat.c
 *
 *  Created on: 17 באוג׳ 2020
 *      Author: גל
 */
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

typedef size_t* Subgroup;

/*
 * A suggestion for Modularity Matrix B struct.
 */
typedef struct _modmat {
	size_t n;	 	/* n=|V| of network */
	spmat A; 		/* Network adjacency matrix in form of spmat */
	size_t *K; 		/* Compact representation of degrees-product matrix (k_i * k_j / M) */
	size_t M; 		/* Total sum of degrees in the network*/
	Subgroup g; 	/* Array of relevant indices of submatrix */

	void (*free)(struct _modmat B, Subgroup g); /*free all resources used by ModMat instance */

	void (*mult)(const struct _modmat B, Subgroup g, const double *v, double *result);
	/*multiply ModMat with vector v and store result.
	 * Can be extended to submatrix B[g].
	 */

} modMat;

void freeModMat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->g);
	free(B);
}

/*allocate new ModMat of size n*/
modMat *allocateModMat(int n){
	modMat *rep=(modMat*)malloc(sizeof(modMat));
	if (rep==NULL)
		exit(1);

	rep->n=n;

	rep->A=spmat_allocate_list(n);
	if (rep->A==NULL){
		free(rep);
		exit(1);
	}

	rep->K=(double*)malloc(n*sizeof(double));
	if (rep->K==NULL){
		free(rep);
		exit(1);
	}

	rep->g=(Subgroup)malloc(n*sizeof(size_t));

	rep->free=freeModMat;

	return rep;
}

/*read bytes from file or exit if failed. used to replace assert functionality*/
void readFromFile(int *dest,unsigned int num,unsigned int size,FILE* file){
	if (fread(dest,num,size,file)!=size)
		exit(FILE_READ_ERROR);
}

/* Load an open input file into previously allocated ModMat B*/

void loadModMatrix(FILE *input, modMat *B){
	unsigned int n, i, currDeg;
	size_t *neighbours, *k, *g;
	readFromFile(n,sizeof(size_t),1,input); /* read n=|V| of network */
	B=allocateModMat(n);
	k=B->K;
	g=B->g;
	/* Populate B with A, K matrices, and compute M */
	for (i=0; i<n; i++){
		readFromFile(currDeg,sizeof(size_t),1,input);
		*k++=currDeg;
		B->M+=currDeg;
		neighbours=(size_t*)malloc(currDeg*sizeof(size_t));
		if (neighbours==NULL)
			exit(MEM_ALLOC_ERROR);
		readFromFile(neighbours,sizeof(size_t),currDeg,input);
		B->A->add_row(B->A,neighbours,i); /*spmat's add_row must be modified to translate indices of V into 0/1 values in A */
		*g++=1;
		free(neighbours);
	}
	rewind(input);
}

/* Multiply submatrix B[g] with vector u;
 * This is an implementation for mult method in ModMat type struct.
 */
void multSubgroupMatrix(modMat B, Subgroup g, double *u, double *res){
	double *resAu, *resKu, sumK=0;
	size_t *j, *i, t, *tempg=g;

	resAu=(double*)malloc(B->n*sizeof(double));
	B->A->mult(B->A,u,g,resAu); /* spmat's mult must be modified to support multiplying B[g] */

	resKu=(double*)malloc(B->n*sizeof(double));
	for (j=B->K; j<B->K+B->n; j++)
		if (*tempg==j-B->K && *tempg<B->n){
			sumK += (*j);
			tempg++;
		}
	/* For original B of the whole network, sumK == B->M */

	tempg=g;
	for (i=B->K; i<B->K+B->n; i++)
		if (*tempg==i-B->K && *tempg<B->n){
			*resKu++ = ((double)((*i)*sumK) / B->M)*(*u++);
			tempg++;
		}
	resKu-=B->n;

	for (t=0; t<B->n; t++)
		*res++=(*resAu++)-(*resKu++);
	free(resAu);
	free(resKu);
}
