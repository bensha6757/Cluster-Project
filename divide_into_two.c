#include "divide_into_two.h"


/*Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s
 * w.r.t to a moved vertex (or INITIAL_Q for initial state).
 * Reuses the vector B[g]_hat * s to enhance computations, if needed.
 */
double get_modularity(modMat *B, vector s, vector Bs, int movedVertex){
	double Q, *p;
	int sgn=0;
	vector tmp;
	if (movedVertex==INITIAL_Q){
		mult_B_hat_g(B, B, s, Bs, false);
		#ifdef DEBUG
		printf("BEGIN: get_modularity_if_multBhat\n");
		#endif
		Q = 0.5 * dot_prod(Bs,s,B->gSize);
	}
	else {
		tmp=(double*)malloc(B->gSize*sizeof(double));
		VERIFY(tmp!=NULL,MEM_ALLOC_ERROR)
		B->get_row(B,movedVertex,tmp);
		*(s+movedVertex)*=-1;
		sgn=*(s+movedVertex);
		for (p=Bs; p<Bs+B->gSize; p++)
			*p=sgn*2*(*(tmp++));
		Q = 0.5 * dot_prod(Bs,s,B->gSize);
		/* Restore s, Bs to initial state */
		*(s+movedVertex)*=-1;
		tmp-=B->gSize;
		for (p=Bs; p<Bs+B->gSize; p++)
			*p+=sgn*2*(*(tmp++));
		free(tmp);
	}
	return Q;
}

/*Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q */
void optimize_division_modified(modMat *B, vector s){
	num i;
	vector max_v=s;/* impInd=0;*/
	double Q_0, Qmax, Qtmp, deltaQ=1,improveTmp=0, improveMax=0;
	double *Bs=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)
	while (deltaQ>0.0){
		Q_0=get_modularity(B ,s, Bs, INITIAL_Q);
		for (i=0; i<B->gSize ; i++){
				Qtmp=get_modularity(B,s,Bs,i)-Q_0;
				if (Qtmp>Qmax){
					Qmax=Qtmp;
					improveTmp+=Qmax;
					max_v=s+i;
				}
				if (improveTmp>improveMax){
					improveMax=improveTmp;
					/*impInd=p-s;*/
				}
		}
		*max_v*=-1;
		deltaQ=improveMax;
	}
	free(Bs);
}

void optimize_division_original(modMat *B, vector s){
	vector p;
	num i, *indices, maxi, *moved, *q, impInd=0;
	double Q_0, *Bs, *score, maxScore=0, deltaQ=1,*improve, improveMax=0;
	
	#ifdef DEBUG
	printf("Begin: Optimize\n");
	#endif
	
	Bs=(double*)malloc(B->gSize*sizeof(double));
	VERIFY (Bs!=NULL,MEM_ALLOC_ERROR)
	score=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(score!=NULL,MEM_ALLOC_ERROR)

	while (deltaQ>0.0){
		moved=(int_vector)calloc(B->gSize,sizeof(num));
		VERIFY(moved!=NULL,MEM_ALLOC_ERROR)
		for (i=0; i<B->gSize ; i++){
			Q_0=get_modularity(B ,s, Bs, INITIAL_Q);
			for (q=moved; q<moved+B->gSize; q++){
				if (!(*q))
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
		for (q=indices+B->gSize; q>indices+impInd; q--)
			*(s+*(q))*=-1;
		deltaQ = q == indices+B->gSize ? 0 : improveMax;
	}
	#ifdef DEBUG
	printf("end1: Optimize\n");
	#endif
	free(Bs);
	free(score);
	free(improve);
	free(indices);
	#ifdef DEBUG
	printf("End: Optimize\n");
	#endif
}

/*
 * Map a real vector to a pre-allocated vector of {-1,1} elements.
 */
void eigen_to_s(modMat *B, vector eigenVec, vector s){
	vector e=eigenVec;
	#ifdef DEBUG
	printf("BEGIN: eigen_to_s\n");
	#endif
	while(e < eigenVec+B->gSize)
		*s++ = IS_POSITIVE(*(e++)) ? 1 : -1;
	#ifdef DEBUG
	printf("SUCCESS: eigen_to_s\n");
	#endif
}

/* Maps a {-1,1} vector s of B's dim. to a partition g1,g2 */  
void map_s_to_groups(modMat *B, vector s, Subgroup *g1, Subgroup *g2,  num *sizeG1, num *sizeG2){
	vector i;
	num v=0;
	*sizeG1=0;
	*sizeG2=0;
	for (i=s; i<s+B->gSize; i++){
		if (IS_POSITIVE(*i))
			v++;
	}
	*sizeG1 = v;
	*sizeG2 = B->gSize - *sizeG1;
	*g1=(int_vector)malloc(*sizeG1*sizeof(num));
	VERIFY(g1!=NULL,MEM_ALLOC_ERROR)
	*g2=(int_vector)malloc(*sizeG2*sizeof(num));
	VERIFY(g2!=NULL,MEM_ALLOC_ERROR)
	for (i=s; i<s+B->gSize; i++){
		v=i-s;
		if (IS_POSITIVE(*i))
			**g1++=v;
		else
			**g2++=v;
	}
}

DIV_RESULT div_into_two(modMat *B,Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	double beta;
	vector u,s;
	modMat *Bg;
	DIV_RESULT ret;
	Bg=create_Sub_Matrix(B, g, sizeG, USE_LINKED);
	u=(vector)malloc(Bg->gSize*sizeof(double));
	VERIFY(u!=NULL,MEM_ALLOC_ERROR)
	leading_eigenpair(B,Bg,u,&beta);
	if (!IS_POSITIVE(beta))
		ret=GROUP_INDIVISIBLE;
	s=(vector)malloc(B->gSize*sizeof(double));
	VERIFY(s!=NULL,MEM_ALLOC_ERROR)
	eigen_to_s(B,u,s);
	if (!IS_POSITIVE(get_modularity(B,s,u,INITIAL_Q)))
		ret=GROUP_INDIVISIBLE;
	else {
		optimize_division_original(B, s);
		map_s_to_groups(B,s, g1, g2, sizeG1, sizeG2);
		ret=GROUP_DIVIDED;
	}
	Bg->free(Bg);
	free(u);
	free(s);
	return ret;
}