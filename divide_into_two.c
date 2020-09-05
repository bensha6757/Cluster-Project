#include "Divide_Into_Two.h"

/** Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
 * 	For enhanced computation, can use an optional pre-computed vector storing the product B*s.
 *  @param given_Bs - boolean; If TRUE, assumes the argument Bs already stores the product B*s.
 * 					  If FALSE, stores the product Bs in a pre-allocated vector Bs.
 */
double get_modularity_init(modMat *B, vector s, vector Bs, boolean given_Bs){
	double Q;
	num gSize=B->gSize;
	if (!given_Bs)
		B->mult(B, s, Bs, NO_SHIFT);
	Q = dot_prod(s,Bs,gSize);

	return Q;
}

/* Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s w.r.t to a moved vertex moved_v.
 * Recieves a pre-computed matrix-vector prodcut B*s, and the vector s itself.
 */
double get_modularity_moved(modMat *B, vector s, vector Bs, num moved_v){
	double Q, *p;
	int sgn=0;
	vector Bi, Bs_mod, Bs_orig;
	num gSize = B->gSize;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: get_modularity_moved for %d\n", moved_v);
	#endif

	Bi=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bi!=NULL,MEM_ALLOC_ERROR)
	B->get_row(B, moved_v, Bi);
	
	Bs_mod=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs_mod!=NULL,MEM_ALLOC_ERROR)

	s[moved_v] *= -1;
	sgn=s[moved_v];

	for (p=Bs_mod, Bs_orig = Bs ; p < Bs_mod + gSize ; p++, Bi++, Bs_orig++)
		*p = *(Bs_orig) + (sgn * 2 * (*Bi));

	Q = dot_prod(Bs_mod,s,gSize);

	/* Restore s to initial state */
	s[moved_v] *= -1;

	free(Bs_mod);
	free(Bi-gSize);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: get_modularity_moved = %f for %d\n",Q, moved_v);
	#endif
	return Q;
}


/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex(modMat *Bg, vector *s, vector Bs, int_vector indices, double *improve) {
	vector p, s_ptr = *s;
	char *m, *moved;
	num i, maxi = 0, gSize = Bg->gSize;
	double Q_0, Q_t, maxScore=0;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: STEP 1 - optimize_division_original\n");
	#endif

	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	moved=(char*)calloc(gSize,sizeof(char));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs, TRUE);
		maxi=0;
		maxScore = -DBL_MAX;
		for (m=moved; m < moved + gSize; m++, p++){
			if (!(*m)) {
				Q_t = get_modularity_moved(Bg, s_ptr, Bs, m-moved) - Q_0;
				if (Q_t > maxScore){
					maxScore = Q_t;
					maxi = m-moved;
				}
			}
		}
		/* lines 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		/* line 19 in Alg. 4 PsCode */
		moved[maxi] = 'm';

		/* lines 14-18 in Alg. 4 PsCode */
		*improve  = (i==0 ? 0 : *(improve-1)) + maxScore;
		improve++;

	}
	free(moved);
	
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division(modMat *Bg, vector *s, vector Bs){
	vector p, s_ptr = *s;
	int_vector indices, j;
	num impInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, *improve, improveMax=0;
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)

	improve=(double*)malloc(gSize*sizeof(double));
	VERIFY(improve!=NULL, MEM_ALLOC_ERROR)

	do {
		lastDeltaQ=deltaQ;
		/* lines 2-20 in Alg. 4 pseudo-code */
		move_maximal_score_vertex(Bg, s, Bs, indices, improve);

		/* line 21 in Alg. 4 pseudo-code */
		improveMax = -DBL_MAX;
		impInd = 0;
		for (p=improve; p < improve + gSize; p++){
			if (*p>improveMax){
				improveMax=*p;
				impInd=p-improve;
			}
		}
		
		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j=indices+gSize-1; j>indices+impInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (impInd == gSize-1) ? 0 : improveMax;

		/* VERIFY(iter < gSize, INFINITE_LOOP_ERROR) */
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && iter++ < gSize);
	#ifndef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division for gSize=%d after %d iter with dQ=%f.\n",gSize, iter, deltaQ);
	#endif
	free(improve);
	free(indices);
}


/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_mod(modMat *Bg, vector *s, vector Bs, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector p, s_ptr = *s;
	char *moved, *m;
	num i, maxi = 0, gSize = Bg->gSize;
	double Q_0, Q_t, maxScore=0, tmpImprove=0;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: STEP 1 - optimize_division_original\n");
	#endif

	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	moved=(char*)calloc(gSize,sizeof(char));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs, TRUE);
		maxi=0;
		maxScore = -DBL_MAX;
		for (m=moved; m < moved + gSize; m++, p++){
			if (!(*m)) {
				Q_t = get_modularity_moved(Bg, s_ptr, Bs, m-moved) - Q_0;
				if (Q_t > maxScore){
					maxScore = Q_t;
					maxi = m-moved;
				}
			}
		}
		/* lines 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		/* line 19 in Alg. 4 PsCode */
		moved[maxi] = 'm';
		
		tmpImprove += maxScore;
		if (tmpImprove > *maxImprove){
			*maxImprove = tmpImprove;
			*maxImpInd = i;
		}
	}
	free(moved);
	
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_mod(modMat *Bg, vector *s, vector Bs){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove;
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)


	do {
		lastDeltaQ=deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex_mod(Bg, s, Bs, indices, &maxImprove, &maxImpInd);

		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j=indices+gSize-1; j>indices+maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		/* VERIFY(iter < gSize, INFINITE_LOOP_ERROR) */
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && iter++ < gSize);
	#ifndef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_mod for gSize=%d after %d iter with dQ=%f.\n",gSize, iter, deltaQ);
	#endif
	free(indices);
}



/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_mod_2(modMat *Bg, vector *s, vector Bs, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	num *moved=NULL;
	num i, j, maxi = 0, gSize = Bg->gSize;
	double Q_0, Q_t, maxScore=0, tmpImprove=0;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: STEP 1 - optimize_division_original\n");
	#endif

	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	moved = allocate_hash_set(gSize);

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs, TRUE);
		maxi=0;
		maxScore = -DBL_MAX;
		for (j=0; j<gSize; j++){
			if (!getFlag(moved, gSize, j)) {
				Q_t = get_modularity_moved(Bg, s_ptr, Bs, j) - Q_0;
				if (Q_t > maxScore){
					maxScore = Q_t;
					maxi = j;
				}
			}
		}
		/* lines 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		/* line 19 in Alg. 4 PsCode */
		setFlag(moved, gSize, maxi);
		
		tmpImprove += maxScore;
		if (tmpImprove > *maxImprove){
			*maxImprove = tmpImprove;
			*maxImpInd = i;
		}
	}
	free(moved);
	
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_mod_2(modMat *Bg, vector *s, vector Bs){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove;
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)


	do {
		lastDeltaQ=deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex_mod_2(Bg, s, Bs, indices, &maxImprove, &maxImpInd);

		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j=indices+gSize-1; j>indices+maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		/* VERIFY(iter < gSize, INFINITE_LOOP_ERROR) */
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && iter++ < gSize);
	#ifndef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_mod for gSize=%d after %d iter with dQ=%f.\n",gSize, iter, deltaQ);
	#endif
	free(indices);
}



/*
 * Map a real vector to a pre-allocated vector of {-1,1} elements.
 */
void eigen_to_s(modMat *Bg, vector eigenVec, vector *s){
	vector e = eigenVec, s_i;
	num gSize = Bg->gSize;
	
	*s=(vector)malloc(gSize*sizeof(double));
	VERIFY(s!=NULL, MEM_ALLOC_ERROR)
	s_i = *s;
	while (e < eigenVec + gSize){
		*(s_i++) = IS_POSITIVE(*(e++)) ? 1 : -1;
	}
	
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

/** Maps a {-1,1} vector s of B's dim. to a partition g1,g2 of g,
 *  where g is some subset of 0,1,...,n, s.t. n=|V| of the original network.
 * */  
void map_s_to_groups(Subgroup g, num sizeG, vector s, Subgroup *g1, Subgroup *g2,  num *sizeG1, num *sizeG2){
	vector s_i;
	num v = 0;
	Subgroup g1_ptr, g2_ptr;
	
	for (s_i=s ; s_i < s + sizeG ; s_i++){
		if (IS_POSITIVE(*s_i))
			v++;
	}
	*sizeG1 = v;
	*sizeG2 = sizeG - v;
	if (v == 0 || v == sizeG){
		return;
	}
	*g1 = (Subgroup)malloc(v * sizeof(num));
	VERIFY(*g1!=NULL,MEM_ALLOC_ERROR)
	*g2 = (Subgroup)malloc((sizeG - v) * sizeof(num));
	VERIFY(*g2!=NULL,MEM_ALLOC_ERROR)

	g1_ptr = *g1;
	g2_ptr = *g2;
	for (s_i = s ; s_i < s + sizeG ; s_i++, g++){
		if (IS_POSITIVE(*s_i)){
			*(g1_ptr++) = *g;
		}
		else {
			*(g2_ptr++) = *g;
		}
	}
	
}


DIV_RESULT divide_into_two(modMat *B, Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	scalar beta;
	vector u, s;
	modMat *Bg;
	DIV_RESULT ret;

	Bg = create_Sub_Matrix(B, g, sizeG);

	leading_eigenpair(Bg, &u, &beta);

	if (!IS_POSITIVE(beta)){
		*sizeG1 = sizeG;
		*sizeG2 = 0;
		return GROUP_INDIVISIBLE;
	}

	eigen_to_s(Bg, u, &s);

	if (!IS_POSITIVE(get_modularity_init(Bg, s, u, FALSE))) /* u will store the product Bs for reuse in optimization (Alg. 4) */
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;

	optimize_division_mod_2(Bg, &s, u);
	Bg->free(Bg);

	map_s_to_groups(g, sizeG, s, g1, g2, sizeG1, sizeG2);

	free(u);
	free(s);
	
	return ret;
}