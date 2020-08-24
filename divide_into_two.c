/*
 * divide_into_two.c
 *
 *  Created on: 24 ����� 2020
 *      Author: ��
 */


#include "io_mem_errors.h"
#include "modMat.h"
#include "divide_into_two.h"
#include "leading_eigenpair.h"
#include <stdio.h>
#include <stdlib.h>

/*Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s */
double get_modularity(modMat *B, vector s, vector Bs, int movedVertex){
	double Q, *p;
	int sgn=0;
	if (movedVertex==INITIAL_Q){
		B->mult(B, B->g, s, Bs);
		Q = 0.5 * dotProd(Bs,s,B->gSize);
	}
	else {
		double *tmp=(double*)malloc(B->gSize*sizeof(double));
		if (tmp==NULL)
			exit(MEM_ALLOC_ERROR);
		B->get_row(B,movedVertex,tmp);
		*(s+movedVertex)*=-1;
		sgn=*(s+movedVertex);
		for (p=Bs; p<Bs+B->gSize; p++)
			*p=sgn*2*(*(tmp++));
		Q = 0.5 * dotProd(Bs,s,B->gSize);
		/* Restore s, Bs to initial state */
		*(s+movedVertex)*=-1;
		tmp-=B->gSize;
		for (p=Bs; p<Bs+B->gSize; p++)
			*p+=sgn*2*(*(tmp++));
		free(tmp);
	}
	return Q;
}


void map_s_to_groups(modMat *B, vector s, size_t *g1, size_t *g2,  size_t *sizeG1, size_t *sizeG2){
	vector i;
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
void optimize_division_modified(modMat *B, vector s){
	vector p;
	size_t i, *max_v=s;/* impInd=0;*/
	double Q_0, Qmax, Qtmp, deltaQ=1,improveTmp=0, improveMax=0;
	double *Bs=(double*)malloc(B->gSize*sizeof(double));
	if (Bs==NULL)
		exit(MEM_ALLOC_ERROR);
	while (deltaQ>0.0){
		Q_0=getModularity(B ,s, Bs, INITIAL_Q);
		for (i=0; i<B->gSize ; i++){
				Qtmp=getModularity(B,s,Bs,i)-Q_0;
				if (Qtmp>Qmax){
					Qmax=Qtmp;
					improveTmp+=Qmax;
					max_v=s+i;
				}
				if (improveTmp>improveMax){
					improveMax=improveTmp;
					/*impInd=p-s;*/
				}
				(*p)*=-1;
		}
		*max_v*=-1;
		deltaQ=improveMax;
	}
	free(Bs);
}

void optimize_division_original(modMat *B, vector s){
	vector p;
	size_t i, *indices, maxi, *moved, impInd=0;
	double Q_0, *score, maxScore=0, deltaQ=1,*improve, improveMax=0;
	double *Bs=(double*)malloc(B->gSize*sizeof(double));
	if (Bs==NULL)
		exit(MEM_ALLOC_ERROR);
	score=(double*)malloc(B->gSize*sizeof(double));
	if (score==NULL){
		free(Bs);
		exit(MEM_ALLOC_ERROR);
	}
	while (deltaQ>0.0){
		moved=(size_t*)calloc(B->gSize,sizeof(size_t));
		if (moved==NULL){
			free(score);
			free(Bs);
			exit(MEM_ALLOC_ERROR);
		}
		for (i=0; i<B->gSize ; i++){
			Q_0=get_modularity(B ,s, Bs, INITIAL_Q);
			for (p=moved; p<moved+B->gSize; p++){
				if (!p)
					*score++=get_modularity(B,s,Bs,i)-Q_0;
			}
			score-=B->gSize;
			for (p=score; p<score+B->gSize; p++){
				if (*p > maxScore)
					maxScore=*p;
					maxi=p-score;
			}
			*(s+maxi)*=-1;
			*indices++  = maxi;
			*improve++ += maxScore;
			*(moved+maxi)=!*(moved+maxi);
		}
		indices-=B->gSize;
		improve-=B->gSize;
		for (p=improve; p<improve+B->gSize; p++){
			if (*p>improveMax){
				improveMax=*p;
				impInd=improve-p;
			}
		}
		for (i=B->gSize-1; i>impInd; i--)
			*(s+*(indices+i))*=-1;
		deltaQ = i == B->gSize-1 ? 0 : improveMax;
	}
	free(Bs);
	free(score);
	free(improve);
	free(indices);
}

/*
 * Map a double's vector to a {-1,1} vector with size_t elements.
 * consider s vector to be a double as well if upcasting fails when applying mult.
 */
void eigen_to_s(modMat *B, vector u, vector s){
	while(u<u+B->gSize)
		*s++ = IS_POSITIVE(*(u++)) ? 1 : -1;
}

int div_into_two(modMat *B,Subgroup g, size_t sizeG, Subgroup **g1, Subgroup **g2, size_t *sizeG1, size_t *sizeG2){
	double beta;
	vector u,s,tmp;
	leading_eigenpair(B,g,u,&beta);
	if (beta<=0.0)
		return GROUP_INDIVISIBLE;
	s=(vector)malloc(B->gSize*sizeof(size_t));
	if (s==NULL)
		exit(MEM_ALLOC_ERROR);
	eigen_to_s(B,u,s);
	tmp=(vector)malloc(B->gSize*sizeof(size_t));
	if (tmp==NULL)
		exit(MEM_ALLOC_ERROR);
	if (get_modularity(B,s,tmp,INITIAL_Q)<=0.0)
		return GROUP_INDIVISIBLE;
	else {
		optimize_division_original(B, s);
		map_s_to_groups(B,s, *g1, *g2, sizeG1, sizeG2);
		free(s);
	}
	return GROUP_DIVIDED;
}

