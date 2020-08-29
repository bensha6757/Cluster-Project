#include "divide_into_two.h"
#define DEBUG

/* Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s
 * w.r.t to a moved vertex (or INITIAL_Q for initial state, before moving any vertex).
 * Reuses the vector B[g]_hat * s to enhance computations, if needed.
 */
double get_modularity(modMat *B, vector s, vector *Bs, int movedVertex){
	double Q, *p;
	int sgn=0;
	vector tmp;
	if (movedVertex==INITIAL_Q){
		mult_B_hat_g(B, B, s, *Bs, false);
		Q = 0.5 * dot_prod(*Bs,s,B->gSize);
	}
	else {
		tmp=(double*)malloc(B->gSize*sizeof(double));
		VERIFY(tmp!=NULL,MEM_ALLOC_ERROR)
		B->get_row(B,movedVertex,tmp);
		*(s+movedVertex)*=-1;
		sgn=*(s+movedVertex);
		for (p=*Bs; p<*Bs+B->gSize; p++)
			*p+=sgn*2*(*(tmp++));
		Q = 0.5 * dot_prod(*Bs,s,B->gSize);
		/* Restore s, Bs to initial state */
		tmp-=B->gSize;
		*(s+movedVertex)*=-1;
		sgn=*(s+movedVertex);
		for (p=*Bs; p<*Bs+B->gSize; p++)
			*p+=sgn*2*(*(tmp++));
		tmp-=B->gSize;
		free(tmp);
	}
	return Q;
}

/*Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q */
void optimize_division_modified(modMat *B, vector *s){
	num i;
	vector max_v=*s;
	double Q_0, Qmax, Qtmp, deltaQ, improveTmp=0, improveMax=0;
	vector Bs=(vector)malloc(B->gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)
	do {
		Q_0=get_modularity(B , *s, &Bs, INITIAL_Q);
		for (i=0; i<B->gSize ; i++){
				Qtmp=get_modularity(B, *s, &Bs, i)-Q_0;
				if (Qtmp>Qmax){
					Qmax=Qtmp;
					improveTmp+=Qmax;
					max_v=*s+i;
				}
				if (improveTmp>improveMax)
					improveMax=improveTmp;
		}
		*max_v*=-1;
		deltaQ=improveMax;
	} while (IS_POSITIVE(deltaQ));
	free(Bs);
}

/* Based on lines 2-20 in Alg. 4 pseudo-code */
void move_maximal_score_vertex(modMat *B, vector *s, int_vector indices, double *improve) {
	vector p, Bs;
	int_vector moved;
	num i, maxi, *m;
	double Q_0, *score, maxScore=0;
	
	Bs=(vector)malloc(B->gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)
	
	/*  A hash set of boolean values, s.t. 
	 *	moved[v]==1 iff vertex v has moved to the opposite group.
	 *  
	 * 	line 2 in Alg. 4 PsCode 
	 */
	moved=(int_vector)calloc(B->gSize,sizeof(num));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	score=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(score!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<B->gSize ; i++){
		/* lines 3-10 in Alg. 4 PsCode */
		Q_0=get_modularity(B ,*s, &Bs, INITIAL_Q);
		p=score;
		for (m=moved; m<moved+B->gSize; m++)
			if (!(*m))
				*p++=get_modularity(B, *s, &Bs, m-moved)-Q_0;
		/* line 11 in Alg. 4 PsCode */
		for (p=score; p<score+B->gSize; p++)
			if (*p > maxScore){
				maxScore=*p;
				maxi=p-score;
			}
		/* lines 12-19 in Alg. 4 PsCode */
		*(*s+maxi)*=-1;
		*indices++  = maxi;
		*improve++ += maxScore;
		*(moved+maxi)=!*(moved+maxi);
	}
	free(moved);
	free(score);
	free(Bs);
	#ifdef DEBUG
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}

/* Based on Alg. 4 pseudo-code */
void optimize_division_original(modMat *B, vector *s){
	vector p;
	int_vector indices, q;
	num impInd=0, iter=0;
	double deltaQ, *improve, improveMax=0;
	#ifdef DEBUG
	printf("BEGIN: optimize_division_original\n");
	#endif
	
	indices=(int_vector)malloc(B->gSize*sizeof(num));
	VERIFY(indices!=NULL,MEM_ALLOC_ERROR)

	improve=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(improve!=NULL,MEM_ALLOC_ERROR)

	do {
		/* lines 2-20 in Alg. 4 pseudo-code */
		move_maximal_score_vertex(B, s, indices, improve);
		/* line 21 in Alg. 4 pseudo-code */
		for (p=improve; p<improve+B->gSize; p++){
			if (*p>improveMax){
				improveMax=*p;
				impInd=p-improve;
			}
		}
		#ifdef DEBUG
		printf("SUCCESS: STEP 2 - optimize_division_original\n");
		#endif
		/* lines 22-25 in Alg. 4 pseudo-code */
		for (q=indices+B->gSize; q>indices+impInd; q--)
			*(*s+*q) *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (impInd == B->gSize-1) ? 0 : improveMax;
		#ifdef DEBUG
		printf("SUCCESS: STEP 3 - optimize_division_original\n");
		#endif
		VERIFY(iter++<B->gSize*B->gSize,INFINITE_LOOP_ERROR)
	} while (IS_POSITIVE(deltaQ));
	free(improve);
	free(indices);
	#ifdef DEBUG
	printf("SUCCESS: optimize_division_original\n");
	#endif
}

/*
 * Map a real vector to a pre-allocated vector of {-1,1} elements.
 */
void eigen_to_s(modMat *B, vector eigenVec, vector *s){
	vector e=eigenVec, p;
	
	#ifdef DEBUG
	printf("BEGIN: eigen_to_s\n");
	#endif
	*s=(vector)malloc(B->gSize*sizeof(double));
	VERIFY(s!=NULL, MEM_ALLOC_ERROR)
	p=*s;
	while (e < eigenVec + B->gSize)
		*p++ = IS_POSITIVE(*(e++)) ? 1 : -1;
	#ifdef DEBUG
	printf("SUCCESS: eigen_to_s\n");
	#endif
}

/* Maps a {-1,1} vector s of B's dim. to a partition g1,g2 */  
void map_s_to_groups(modMat *B, vector s, Subgroup *g1, Subgroup *g2,  num *sizeG1, num *sizeG2){
	vector i;
	num v=0;
	#ifdef DEBUG
	printf("BEGIN: map_s_to_groups\n");
	#endif
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
			(*(*g1)++)=v;
		else
			(*(*g2)++)=v;
	}
	#ifdef DEBUG
	printf("SUCCESS: map_s_to_groups\n");
	#endif
}

DIV_RESULT div_into_two(modMat *B,Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	scalar beta;
	vector u, s;
	modMat *Bg;
	DIV_RESULT ret;

	#ifdef DEBUG
	printf("BEGIN: div_into_two\n");
	#endif

	Bg=create_Sub_Matrix(B, g, sizeG, USE_LINKED);
	leading_eigenpair(B, Bg, &u, &beta);
	if (!IS_POSITIVE(beta))
		ret=GROUP_INDIVISIBLE;
	eigen_to_s(Bg, u, &s);
	if (!IS_POSITIVE(get_modularity(Bg, s, &u, INITIAL_Q)))
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;
	optimize_division_original(Bg, &s);
	map_s_to_groups(Bg, s, g1, g2, sizeG1, sizeG2);
	Bg->free(Bg);
	free(u);
	free(s);
	#ifdef DEBUG
	printf("SUCCESS: div_into_two\n");
	#endif
	return ret;
}