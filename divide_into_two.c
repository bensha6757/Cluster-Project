#include "divide_into_two.h"

/*  Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
	As a secondary result, assigns the product vector B*s to Bs argument.
 */

double get_modularity_init(modMat *Bg, vector s, vector Bs){
	double Q;
	num gSize=Bg->gSize;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: get_modularity_init\n");
	#endif

	Bg->mult(Bg, s, Bs, NO_SHIFT);
	Q = dot_prod(Bs,s,gSize);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: get_modularity_init = %f\n",Q);
	#endif
	return Q;
}

/* Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s w.r.t to a moved vertex moved_v.
 * Recieves a pre-computed matrix-vector prodcut B*s, and the vector s itself.
 */
double get_modularity_moved(modMat *Bg, vector s, vector Bs, num moved_v){
	double Q, *p;
	int sgn=0;
	vector Bi, Bs_mod, Bs_orig;
	num gSize = Bg->gSize;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: get_modularity_moved for %d\n", moved_v);
	#endif

	Bi=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bi!=NULL,MEM_ALLOC_ERROR)
	Bg->get_row(Bg,moved_v,Bi);
	
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
void move_maximal_score_vertex(modMat *Bg, vector *s, int_vector indices, double *improve) {
	vector p, Bs, s_ptr = *s;
	int_vector moved;
	num i, maxi = 0, *m, gSize = Bg->gSize;
	double Q_0, Q_t, *score, maxScore=0;

	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: STEP 1 - optimize_division_original\n");
	#endif

	/*  line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 *	moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	moved=(int_vector)calloc(gSize,sizeof(num));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	score=(double*)malloc(gSize*sizeof(double));
	VERIFY(score!=NULL,MEM_ALLOC_ERROR)

	Bs=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<gSize; i++){
		/* lines 3-10 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs);
		for (p=score, m=moved; m < moved + gSize; m++, p++){
			if (!(*m)) {
				Q_t = get_modularity_moved(Bg, s_ptr, Bs, m-moved) - Q_0;
				*p=Q_t;
				/*
				s_ptr[m-moved]*=-1;
				*p = get_modularity_init(Bg, s_ptr, Bs) - Q_0;
				s_ptr[m-moved]*=-1;
				*/
			}
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
	
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - optimize_division_original\n");
	#endif
}

/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_modified(modMat *Bg, vector *s, int_vector indices, double *improve) {
	vector p, Bs, s_ptr = *s;
	int_vector moved;
	num i, maxi = 0, *m, gSize = Bg->gSize;
	double Q_0, Q_t, maxScore=0;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: STEP 1 - optimize_division_original\n");
	#endif
	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	moved=(int_vector)calloc(gSize,sizeof(num));
	VERIFY(moved!=NULL,MEM_ALLOC_ERROR)

	Bs=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs!=NULL,MEM_ALLOC_ERROR)

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs);
		maxi=0;
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
		/* lines 14-18 in Alg. 4 PsCode */
		*improve  = (i==0 ? 0 : *(improve-1)) + maxScore;
		improve++;
		/* line 19 in Alg. 4 PsCode */
		moved[maxi] = 1;
		maxScore=*(improve-i);
	}
	free(Bs);
	free(moved);
	
	#ifdef DEBUG_DIV_TWO
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
	
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: optimize_division_original!!!\n");
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
		#ifdef DEBUG_DIV_TWO
		printf("SUCCESS: STEP 2 - optimize_division_original\n");
		#endif
		/* lines 22-25 in Alg. 4 pseudo-code */
		for (q=indices+gSize-1; q>indices+impInd; q--)
			s_ptr[*q] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (impInd == gSize-1) ? 0 : improveMax;
		#ifdef DEBUG_DIV_TWO
		printf("SUCCESS: STEP 3 - optimize_division_original, deltaQ=%f\n",deltaQ);
		#endif
		iter++;
		VERIFY(iter++ < gSize * gSize, INFINITE_LOOP_ERROR)
	} while (IS_POSITIVE(deltaQ));
	free(improve);
	free(indices);
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_original\n");
	#endif
}

/*
 * Map a real vector to a pre-allocated vector of {-1,1} elements.
 */
void eigen_to_s(modMat *Bg, vector eigenVec, vector *s){
	vector e = eigenVec, p;
	num gSize = Bg->gSize;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: eigen_to_s\n");
	#endif
	*s=(vector)malloc(gSize*sizeof(double));
	VERIFY(s!=NULL, MEM_ALLOC_ERROR)
	p = *s;
	while (e < eigenVec + gSize){
		*(p++) = IS_POSITIVE(*(e++)) ? 1 : -1;
	}
	#ifdef DEBUG_DIV_TWO
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

/** Maps a {-1,1} vector s of B's dim. to a partition g1,g2 of g,
 *  where g is some subset of 0,1,...,n, s.t. n=|V| of the original network.
 * */  
void map_s_to_groups(Subgroup g, num sizeG, vector s, Subgroup *g1, Subgroup *g2,  num *sizeG1, num *sizeG2){
	vector i;
	num v = 0;
	Subgroup g1_ptr, g2_ptr;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: map_s_to_groups\n");
	#endif
	for (i=s ; i < s + sizeG ; i++){
		if (IS_POSITIVE(*i))
			v++;
	}
	*sizeG1 = v;
	*sizeG2 = sizeG - v;
	*g1 = (Subgroup)malloc(v * sizeof(num));
	VERIFY(*g1!=NULL,MEM_ALLOC_ERROR)
	*g2 = (Subgroup)malloc((sizeG - v) * sizeof(num));
	VERIFY(*g2!=NULL,MEM_ALLOC_ERROR)

	g1_ptr = *g1;
	g2_ptr = *g2;
	for (i=s ; i < s + sizeG ; i++, g++){
		if (IS_POSITIVE(*i)){
			*(g1_ptr++) = *g;
		}
		else{
			*(g2_ptr++) = *g;
		}
	}
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: map_s_to_groups\n");
	#endif
}


DIV_RESULT div_into_two(modMat *B, Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	scalar beta;
	vector u, s;
	modMat *Bg;
	DIV_RESULT ret;

	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: div_into_two\n");
	printG2(g,sizeG);
	#endif

	Bg = create_Sub_Matrix(B, g, sizeG);

	leading_eigenpair(Bg, &u, &beta);

	if (!IS_POSITIVE(beta)){
		*sizeG1 = sizeG;
		*g1=(Subgroup)malloc(sizeof(num)*sizeG);
		VERIFY(*g1!=NULL,MEM_ALLOC_ERROR)
		memcpy(*g1, g, sizeof(num)*sizeG);
		*sizeG2 = 0;
		return GROUP_INDIVISIBLE;
	}

	eigen_to_s(Bg, u, &s);

	if (!IS_POSITIVE(get_modularity_init(Bg, s, u))) /* u is not needed anymore at this stage */
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;
	
	#ifdef DEBUG_DIV_TWO
	printS(s,Bg->gSize);
	#endif

	optimize_division_original(Bg, &s);
	Bg->free(Bg);
	#ifdef DEBUG_DIV_TWO
	printS(s,Bg->gSize);
	#endif

	map_s_to_groups(g, sizeG, s, g1, g2, sizeG1, sizeG2);

	free(u);
	free(s);
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: div_into_two\n");
	#endif
	return ret;
}