/*
 * algorithm.c
 *
 *  Created on: 14 באוג׳ 2020
 *      Author: גל
 */

#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>

#define IS_POSITIVE(X) ((X) > 0.00001)

typedef enum message_t {
	GROUP_INDIVISIBLE,
	GROUP_DIVIDED
} MESSAGE;

typedef size_t* s_vector;

typedef double* vector;
typedef unsigned int boolean;

double dotProd(vector v, vector u, size_t d){
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



void powerIteration(modMat *B, modMat *Bg, vector v, vector result){
	vector p;
	double nrm;
	Bg->mult(B, Bg,v,result);
	nrm=norm(result,Bg->gSize);
	for (p=result; p<result+Bg->gSize; p++)
		*p/=nrm;
}

void setRandVector(vector v, size_t n){
	unsigned int i;
	for (i=0; i<n; i++, v++)
		*v = (double)rand();
}

/*	Approximate dominant eigen value of matrix B
 *  with last vector v returned from power iterations.*/
double approxDomEigenVal(modMat *B, vector v){
	return B->mult(B,v) / dotProd(v,v,B->gSize);
}

boolean isWithin(vector a, vector b, size_t d){
	unsigned int i;
	for (i=0; i<d; i++, a++, b++){
		if (IS_POSITIVE(fabs(*a - *b)))
			return 0;
	}
	return 1;
}


void leadingEigenPair(modMat *B, modMat *Bg, vector leadEigenVec, double* leadEigenVal){
	size_t iter=0;
	vector bprev=(vector)malloc(Bg->gSize*sizeof(double));
	if (bprev==NULL)
		exit(MEM_ALLOC_ERROR);
	vector bnext=(vector)malloc(Bg->gSize*sizeof(double));
	if (bnext==NULL){
		free(bprev);
		exit(MEM_ALLOC_ERROR);
	}
	setRandVector(bprev, Bg->gSize);
	powerIteration(B,Bg,bprev,bnext);
	iter++;
	while (!isWithin(bprev,bnext,Bg->gSize)){
		memcpy(bprev,bnext,Bg->gSize*sizeof(double*));
		powerIteration(B,Bg,bprev,bnext);
		iter++;
	}
	#ifdef PERFORMANCE_ITER
		printf("# of power iterations: %d\n", (int)iter);
	#endif
	leadEigenVec=bnext;
	*leadEigenVal = approxDomEigenVal(Bg, bnext) - getOneNorm(Bg);
	free(bnext);
	free(bprev);
}

/*Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s */
double getModularity(modMat *B, Subgroup g, s_vector s){
	double Q;
	double *Bs=(double*)malloc(B->gSize*sizeof(double));
	B->mult(B, B->g, s, Bs);
	Q = 0.5 * dotProd(Bs,s,B->gSize);
	free(Bs);
	return Q;
}

void mapSToGroups(modMat *B, size_t *s, size_t *g1, size_t *g2,  size_t *sizeG1, size_t *sizeG2){
	size_t *i;
	size_t v=0;
	*sizeG1=0;
	*sizeG2=0;
	for (i=s; i<s+B->gSize; i++){
		if (IS_POSITIVE(*i))
			v++;
	}
	*sizeG1 = v;
	*sizeG2 = B->gSize - *sizeG1;
	g1=(size_t*)malloc(*sizeG1*sizeof(size_t));
	if (g1==NULL)
		exit(MEM_ALLOC_ERROR);
	g2=(size_t*)malloc(*sizeG2*sizeof(size_t));
	if (g2==NULL){
		free(g1);
		exit(MEM_ALLOC_ERROR);
	}
	for (i=s; i<s+B->gSize; i++){
		v=i-s;
		if (IS_POSITIVE(*i))
			*g1++=v;
		else
			*g2++=v;
	}

}

/*Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q */
void optimizeDivision(modMat *B, s_vector s){
	s_vector entry=s;
	size_t *max_v=0, *indices;
	double Qmax, Qtmp deltaQ=1;
	while (deltaQ>0.0){
		for (;entry<s+B->gSize ; entry++){
				(*entry)*=-1; /*Temporarely move vertex s[i] to the other group */
				Qnext=getModularity(B,B->g,s);
				if (Qnext>=Qprev){
					Qprev=Qnext;
					max_v=*entry;
				}
				(*entry)*=-1;
			}
			*max_v*=-1;
	}
	Qprev=getModularity(B,B->g,s);

}

void eigenToS(modMat *B, vector u, s_vector s){
	while(u<u+B->gSize)
		*s++ = IS_POSITIVE(*u++) ? 1 : -1;
}

int div_into_two(modMat *B,Subgroup g, size_t sizeG, Subgroup **g1, Subgroup **g2, size_t *sizeG1, size_t *sizeG2){
	double beta;
	vector u;
	s_vector s;
	leadingEigenPair(B,g,u,&beta);
	if (beta<=0.0)
		return GROUP_INDIVISIBLE;
	s=(s_vector)malloc(B->gSize*sizeof(size_t));
	if (s==NULL)
		exit(MEM_ALLOC_ERROR);
	eigenToS(B,u,s);
	if (getModularity(B,g,s)<=0.0)
		return GROUP_INDIVISIBLE;
	else {
		optimizeDivision(B, s);
		mapSToGroups(B,s, *g1, *g2, sizeG1, sizeG2);
		free(s);
	}
	return GROUP_DIVIDED;
}
