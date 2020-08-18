/*
 * algorithm.c
 *
 *  Created on: 14 באוג׳ 2020
 *      Author: גל
 */

#include "modMat.c"
#define G1 1
#define G2 -1
typedef enum message_t {
	GROUP_INDIVISABLE,
	GROUP_DIVIDED
} MESSAGE;

typedef int* s_vector;


void powerIteration(vector ret, spmat *mat, vector v, dimension d){
	unsigned int i;
	double nrm;
	mat->mult(mat,v,ret);
	nrm=norm(ret,d);
	for (i=0; i<d; i++, ret++)
		*ret/=nrm;
	#ifdef DEBUG
		printf("main: Iteration done\n");
	#endif
}

/*Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s */
double getModularity(ModMat *B,Subgroup g, s_vector s){
	double Q;
	int *res1=(*int)malloc(B->n*sizeof(int));
	double *res2=(*double)malloc(B->n*sizeof(double));
	B->A->mult(B->A, g, s, res1);
	B->K->mult(B->K, g, s, res2); /* consider M */
	Q = 0.5 *(dotProd(res1,s,B->n)-dotProd(res2,s,B->n));
	free(res1);
	free(res2);
	return Q;
}

/*Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q */
void optimizeDivision(ModMat *B, s_vector s){
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

