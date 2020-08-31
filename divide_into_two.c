#include "divide_into_two.h"
#define DEBUG

/*  Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
	As a secondary result, assigns the product vector B*s to Bs argument.
 */

double get_modularity_init(modMat *Bg, vector s, vector *Bs){
	double Q;
	num gSize=Bg->gSize;
	vector Bs_tmp;
	/*printf("BEGIN: get_modularity_init\n");*/
	Bs_tmp=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs_tmp!=NULL,MEM_ALLOC_ERROR)

	mult_B_hat_g(Bg, s, Bs_tmp, NO_SHIFT);
	Q = 0.5 * dot_prod(Bs_tmp,s,gSize);

	memcpy(*Bs,Bs_tmp,gSize*sizeof(double));
	free(Bs_tmp);
	/*printf("SUCCESS: get_modularity_init = %f\n",Q);*/
	return Q;
}

/* Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s w.r.t to a moved vertex moved_v.
 * Recieves a pre-computed matrix-vector prodcut B*s, and the vector s itself.
 */
double get_modularity_moved(modMat *Bg, vector s, vector Bs, num moved_v){
	double Q, *p;
	int sgn=0;
	vector Bi_tmp, Bs_tmp, s_tmp, e_i;
	num gSize = Bg->gSize;
	/*printf("BEGIN: get_modularity_moved for %d\n", moved_v);*/
	Bi_tmp=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bi_tmp!=NULL,MEM_ALLOC_ERROR)
	Bs_tmp=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs_tmp!=NULL,MEM_ALLOC_ERROR)
	s_tmp=(vector)malloc(gSize*sizeof(double));
	VERIFY(s_tmp!=NULL,MEM_ALLOC_ERROR)

	memcpy(s_tmp, s, gSize*sizeof(double));
	memcpy(Bs_tmp, Bs, gSize*sizeof(double));

	/*Bg->get_row(Bg,moved_v,Bi_tmp);*/
	get_basis_unit_vec(&e_i,moved_v,gSize);
	mult_B_hat_g(Bg, e_i, Bi_tmp, NO_SHIFT);
	free(e_i);
	
	*(s_tmp+moved_v) *= -1;
	sgn=*(s_tmp+moved_v);

	for (p=Bs_tmp ; p < Bs_tmp + gSize ; p++)
		*p += sgn * 2 * (*(Bi_tmp++));
	Q = 0.5 * dot_prod(Bs_tmp,s_tmp,gSize);

	free(Bi_tmp-gSize);
	free(Bs_tmp);
	free(s_tmp);
	
	/*printf("SUCCESS: get_modularity_moved = %f for %d\n",Q, moved_v);*/
	return Q;
}

/*
void optimize_division_modified(modMat *B, modMat *Bg, vector *s){
	num i;
	vector max_v=*s, s_ptr = *s;
	double Q_0, Qmax, Qtmp, deltaQ, improveTmp=0, improveMax=0;
	vector Bs=(vector)malloc(Bg->gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)
	do {
		Q_0=get_modularity_init(B, Bg, s_ptr, &Bs);
		for (i=0; i<Bg->gSize ; i++){
				Qtmp=get_modularity_moved(B, Bg, s_ptr, Bs, i)-Q_0;
				if (Qtmp>Qmax){
					Qmax = Qtmp;
					improveTmp += Qmax;
					max_v = s_ptr + i;
				}
				if (improveTmp>improveMax)
					improveMax=improveTmp;
		}
		*max_v*=-1;
		deltaQ=improveMax;
	} while (IS_POSITIVE(deltaQ));
	free(Bs);
}
*/

/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex(modMat *Bg, vector *s, int_vector indices, double *improve) {
	vector p, Bs, s_ptr = *s;
	int_vector moved;
	num i, maxi = 0, *m, gSize = Bg->gSize;
	double Q_0, *score, maxScore=0;
	
	/*  A hash set of boolean values, s.t. 
	 *	moved[v]==1 iff vertex v has moved to the opposite group.
	 *  
	 * 	line 2 in Alg. 4 PsCode 
	 */
	moved=(int_vector)calloc(gSize,sizeof(num));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	score=(double*)malloc(gSize*sizeof(double));
	VERIFY(score!=NULL,MEM_ALLOC_ERROR)

	Bs=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<gSize; i++){
		/* lines 3-10 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, &Bs);
		p = score;
		for (m=moved; m < moved + gSize; m++, p++){
			if (!(*m)) 
				*p = get_modularity_moved(Bg, s_ptr, Bs, m-moved) - Q_0;
		}
		/* line 11 in Alg. 4 PsCode */
		maxi = 0;
		maxScore = *score;
		m=moved;
		for (p = score; p < score + gSize; p++, m++){
			if (*p > maxScore && !(*m)){
				maxScore = *p;
				maxi = p - score;
			}
		}
		/* lines 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		/* lines 14-18 in Alg. 4 PsCode */
		*improve  = (i==0 ? 0 : *(improve-1)) + maxScore;
		improve++;
		/* line 19 in Alg. 4 PsCode */
		moved[maxi] = 1;
	}
	free(Bs);
	free(score);
	free(moved);
	
	#ifdef DEBUG
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}

/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_original(modMat *Bg, vector *s){
	vector p, s_ptr = *s;
	int_vector indices, q;
	num impInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ, *improve, improveMax=0;
	#ifdef DEBUG
	printf("BEGIN: optimize_division_original\n");
	#endif
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL,MEM_ALLOC_ERROR)

	improve=(double*)malloc(gSize*sizeof(double));
	VERIFY(improve!=NULL,MEM_ALLOC_ERROR)

	do {
		/* lines 2-20 in Alg. 4 pseudo-code */
		move_maximal_score_vertex(Bg, s, indices, improve);
		/* line 21 in Alg. 4 pseudo-code */
		improveMax = *improve;
		impInd = 0;
		for (p=improve; p < improve + gSize; p++){
			if (*p>improveMax){
				improveMax=*p;
				impInd=p-improve;
			}
		}
		#ifdef DEBUG
		printf("SUCCESS: STEP 2 - optimize_division_original\n");
		#endif
		/* lines 22-25 in Alg. 4 pseudo-code */
		for (q=indices+gSize-1; q>indices+impInd; q--)
			s_ptr[*q] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (impInd == gSize-1) ? 0 : improveMax;
		#ifdef DEBUG
		printf("SUCCESS: STEP 3 - optimize_division_original\n");
		#endif
		VERIFY(iter++ < gSize * gSize,INFINITE_LOOP_ERROR)
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
void eigen_to_s(modMat *Bg, vector eigenVec, vector *s){
	vector e = eigenVec, p;
	num gSize = Bg->gSize;
	#ifdef DEBUG
	printf("BEGIN: eigen_to_s\n");
	#endif
	*s=(vector)malloc(gSize*sizeof(double));
	VERIFY(s!=NULL, MEM_ALLOC_ERROR)
	p = *s;
	while (e < eigenVec + gSize){
		*(p++) = IS_POSITIVE(*(e++)) ? 1 : -1;
	}
	#ifdef DEBUG
	printf("SUCCESS: eigen_to_s\n");
	#endif
}

void printG2(Subgroup g, int n){
	Subgroup p;
	for (p = g; p < g + n ; p++){
		printf("%d, ", *p);
	}
	printf("\n");
}

void printS(vector s, int n){
	vector p;
	for (p = s; p < s + n ; p++){
		printf("%d, ", (IS_POSITIVE(*p) ? 1 : -1));
	}
	printf("\n");
}

/* Maps a {-1,1} vector s of B's dim. to a partition g1,g2 */  
void map_s_to_groups(modMat *Bg, vector s, Subgroup *g1, Subgroup *g2,  num *sizeG1, num *sizeG2){
	vector i;
	num v = 0, gSize = Bg->gSize;
	Subgroup g = Bg->g, g1_ptr, g2_ptr;
	#ifdef DEBUG
	printf("BEGIN: map_s_to_groups\n");
	#endif
	for (i=s ; i < s + gSize ; i++){
		if (IS_POSITIVE(*i))
			v++;
	}
	*sizeG1 = v;
	*sizeG2 = gSize - v;
	*g1 = (Subgroup)malloc(v * sizeof(num));
	VERIFY(*g1!=NULL,MEM_ALLOC_ERROR)
	*g2 = (Subgroup)malloc((gSize - v) * sizeof(num));
	VERIFY(*g2!=NULL,MEM_ALLOC_ERROR)

	g1_ptr = *g1;
	g2_ptr = *g2;
	for (i=s ; i < s + gSize ; i++, g++){
		if (IS_POSITIVE(*i)){
			*(g1_ptr++) = *g;
		}
		else{
			*(g2_ptr++) = *g;
		}
	}
	#ifdef DEBUG
	printf("SUCCESS: map_s_to_groups\n");
	#endif
}


DIV_RESULT div_into_two(modMat *B, Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	scalar beta;
	vector u, s;
	modMat *Bg;
	DIV_RESULT ret;
	num gSize;

	#ifdef DEBUG
	printf("BEGIN: div_into_two\n");
	#endif
	printG2(g,sizeG);
	Bg = create_Sub_Matrix(B, g, sizeG, USE_LINKED);
	gSize = Bg->gSize;
	printf("sizeG = %d\n", gSize);
	leading_eigenpair(Bg, &u, &beta);
	if (!IS_POSITIVE(beta))
		ret=GROUP_INDIVISIBLE;
	eigen_to_s(Bg, u, &s);
	if (!IS_POSITIVE(get_modularity_init(Bg, s, &u))) /* u is not needed anymore at this stage */
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;
	printS(s,gSize);
	optimize_division_original(Bg, &s);
	printS(s,gSize);
	map_s_to_groups(Bg, s, g1, g2, sizeG1, sizeG2);
	Bg->free(Bg);
	free(u);
	free(s);
	#ifdef DEBUG
	printf("SUCCESS: div_into_two\n");
	#endif
	return ret;
}