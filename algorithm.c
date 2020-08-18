/*
 * algorithm.c
 *
 *  Created on: 14 באוג׳ 2020
 *      Author: גל
 */

#include "modMat.h"
#define G1 1
#define G2 -1
typedef enum message_t {
	GROUP_INDIVISABLE,
	GROUP_DIVIDED
} MESSAGE;

typedef int* s_vector;
typedef double* vector;

void powerIteration(vector ret, modMat *B, vector v){
	unsigned int i;
	double nrm;
	B->mult(B,B->g,v,ret);
	nrm=norm(ret,B->n);
	for (i=0; i<B->n; i++, ret++)
		*ret/=nrm;
}

/*Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s */
double getModularity(modMat *B, Subgroup g, s_vector s){
	double Q;
	double *Bs=(double*)malloc(B->n*sizeof(double));
	B->mult(B, B->g, s, Bs);
	Q = 0.5 * dotProd(Bs,s,B->n);
	free(Bs);
	return Q;
}

/*Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q */
void optimizeDivision(modMat *B, s_vector s){
	int *entry=s;
	double Qprev,Qnext;
	Qprev=B->getModularity(B,B->g,s);
	for (;entry < s+B->n ; entry++){
		(*entry)*=-1; /*Temporarely move vertex s[i] to the other group */
		Qnext=B->getModularity(B,B->g,s);
		if (Qnext>=Qprev)
			Qprev=Qnext;
		else
			(*entry)*=-1;
		/* If Modularity has not decreased, vertex s[i] will be moved to the other group.
		 * Otherwise, s[i] stays in its original group.
		 */
	}
}

