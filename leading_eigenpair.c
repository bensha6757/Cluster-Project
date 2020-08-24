/*
 * leading_eigenpair.c
 *
 *  Created on: 24 ����� 2020
 *      Author: ��
 */

/*
 * algorithm.c
 *
 *  Created on: 14 ����� 2020
 *      Author: ��
 */

#include "io_mem_errors.h"
#include "modMat.h"
#include "leading_eigenpair.h"
#include <stdio.h>
#include <stdlib.h>

double dot_prod(vector v, vector u, size_t d){
	unsigned int i;
	double acc=0;
	for (i=0; i<d; i++, v++, u++){
		acc += (*v)*(*u);
	}
	return acc;
}

double norm(vector v, size_t d){
	return sqrt(dotProd(v,v,d));
}

void set_rand_vector(vector v, size_t n){
	unsigned int i;
	for (i=0; i<n; i++, v++)
		*v = (double)rand();
}

/*	Approximate dominant eigen value of matrix B
 *  with last vector v returned from power iterations.*/
double approx_dom_eigen_val(modMat *B, vector bprev, vector bnext){
	return dot_prod(bnext,bprev,B->gSize) / dot_prod(bprev,bprev,B->gSize);
}

boolean is_within(vector a, vector b, size_t d){
	size_t i;
	for (i=0; i<d; i++, a++, b++){
		if (IS_POSITIVE(fabs(*a - *b)))
			return 0;
	}
	return 1;
}

/*
 * Apply a single iteration to a matrix B and a vector v according to the power method.
 * Store result in result vector.
 */
void power_iteration(modMat *B, modMat *Bg, vector v, vector result){
	vector p;
	double nrm;
	Bg->mult(B, Bg, v, result);
	nrm=norm(result,Bg->gSize);
	for (p=result; p<result+Bg->gSize; p++)
		*p/=nrm;
}

/*
 * Compute leading eigen pair of modularity Matrix B_hat[g]
 */
void leading_eigenpair(modMat *B, modMat *Bg, vector leadEigenVec, double* leadEigenVal){
	size_t iter=0;
	vector bprev=(vector)malloc(Bg->gSize*sizeof(double));
	if (bprev==NULL)
		exit(MEM_ALLOC_ERROR);
	vector bnext=(vector)malloc(Bg->gSize*sizeof(double));
	if (bnext==NULL){
		free(bprev);
		exit(MEM_ALLOC_ERROR);
	}
	set_rand_vector(bprev, Bg->gSize);
	power_iteration(B,Bg,bprev,bnext);
	iter++;
	while (!is_within(bprev,bnext,Bg->gSize)){
		memcpy(bprev,bnext,Bg->gSize*sizeof(double*));
		power_iteration(B,Bg,bprev,bnext);
		iter++;
		if (iter>1000*Bg->gSize)
			exit(INFINITE_LOOP_ERROR);
	}
	#ifdef PERFORMANCE_ITER
		printf("# of power iterations: %d\n", (int)iter);
	#endif
	leadEigenVec=bnext;
	*leadEigenVal = approx_dom_eigen_val(Bg,bprev,bnext) - get_1_norm(Bg);
	free(bnext);
	free(bprev);
}
