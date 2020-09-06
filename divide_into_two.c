#include "Divide_Into_Two.h"

/** Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
 * 	For enhanced computation, can use an optional pre-computed vector storing the product B*s.
 *  @param given_Bs - boolean; If TRUE, assumes the argument Bs already stores the product B*s.
 * 					  If FALSE, stores the product Bs in a pre-allocated vector Bs.
 */
double get_modularity_init(modMat *B, vector s, vector Bs, boolean given_Bs){
	double Q;
	num gSize=B->gSize;
	/*clock_t start, end;*/


	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: get_modularity_init\n");
	#endif

	if (!given_Bs)
		B->mult(B, s, Bs, NO_SHIFT);
	Q = dot_prod(s,Bs,gSize);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: get_modularity_init\n");
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif

	return Q;
}

void modify_Bs(modMat *B, vector Bs, int i, int sgn){
	num gSize=B->gSize;
	vector Bi, p;

	Bi=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bi!=NULL,MEM_ALLOC_ERROR)
	B->get_row(B, i, Bi);

	for (p=Bs ; p < Bs + gSize ; p++, Bi++)
		*p += (sgn * 2 * (*Bi));
}


/* Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s w.r.t to a moved vertex moved_v.
 * Recieves a pre-computed matrix-vector prodcut B*s, and the vector s itself.
 */
double get_modularity_moved(modMat *B, vector s, vector Bs, num moved_v){
	double Q;
	vector Bs_copy;
	num gSize = B->gSize;
	#ifdef DEBUG_DIV_TWO
	printf("BEGIN: get_modularity_moved for %d\n", moved_v);
	#endif

	Bs_copy=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bs_copy!=NULL,MEM_ALLOC_ERROR)
	memcpy(Bs_copy, Bs, gSize*sizeof(double));

	s[moved_v] *= -1;

	modify_Bs(B, Bs_copy, moved_v, (int)s[moved_v]);
	Q = dot_prod(s, Bs_copy, gSize);

	/* Restore s to initial state */
	s[moved_v] *= -1;

	free(Bs_copy);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: get_modularity_moved = %f for %d\n",Q, moved_v);
	#endif
	return Q;
}

double BEN_get_modularity_moved(modMat *Bg, vector s, vector Bs, num moved_v){
	double Q = 0;
	int sgn = 0;
	vector Bi, Bs_orig, s_i;
	num gSize = Bg->gSize;
	/*clock_t start, end;*/

	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: BEN_get_modularity_moved\n");
	#endif

	Bi=(vector)malloc(gSize*sizeof(double));
	VERIFY(Bi!=NULL,MEM_ALLOC_ERROR)
	Bg->get_row(Bg,moved_v,Bi);

	s[moved_v] *= -1;
	sgn = s[moved_v];

	for (s_i = s, Bs_orig = Bs ; Bs_orig < Bs + gSize ; s_i++, Bi++, Bs_orig++){
		Q += ((*Bs_orig + (sgn * 2 * (*Bi))) * (*s_i));
	}

	/* Restore s to initial state */
	s[moved_v] *= -1;

	free(Bi-gSize);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: BEN_get_modularity_moved\n");
	end = clock();
	printf("Execution took %f seconds\n ", ((double)(end-start) / CLOCKS_PER_SEC));
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
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && ++iter < gSize);
	#ifndef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_mod for gSize=%d after %d iter with dQ=%f.\n",gSize, iter, deltaQ);
	#endif
	free(indices);
}



/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_mod_2(modMat *Bg, vector *s, vector Bs, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	long_num *moved=NULL;
	num i, k, maxi = 0, gSize = Bg->gSize;
	double Q_0, Q_t, maxScore = -DBL_MAX, tmpImprove=0;
	/*clock_t start, end;*/

	
	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: STEP 1 - move_maximal_score_vertex_mod_2\n");
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
		for (k=0; k<gSize; k++){
			if (!getFlag(moved, gSize, k)) {
				Q_t = BEN_get_modularity_moved(Bg, s_ptr, Bs, k) - Q_0;
				if (Q_t > maxScore){
					maxScore = Q_t;
					maxi = k;
				}
			}
		}
		/* line 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		modify_Bs(Bg, Bs, maxi , s_ptr[maxi]);

		/* line 14-18, 21 in Alg. 4 PsCode */
		tmpImprove += maxScore;
		if (tmpImprove > *maxImprove){
			*maxImprove = tmpImprove;
			*maxImpInd = i;
		}
		/* line 19 in Alg. 4 PsCode */
		setFlag(moved, gSize, maxi);
	}
	free(moved);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - move_maximal_score_vertex_mod_2\n");
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_mod_2(modMat *Bg, vector *s, vector Bs){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove=-DBL_MAX;
	/* clock_t start, end; */
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)

	
	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: optimize_division_original_mod_2\n");
	#endif
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
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && ++iter < gSize);
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_original_mod_2\n");
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif
	free(indices);
}

/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_mod_3(modMat *Bg, vector *s, vector Bs, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	long_num *unmoved=NULL;
	num i, maxi = 0, gSize = Bg->gSize;
	int k=0;
	double Q_0, Q_t, maxScore = -DBL_MAX, tmpImprove=0;
	/*clock_t start, end;*/

	
	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: STEP 1 - move_maximal_score_vertex_mod_2\n");
	#endif

	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * unmoved[v]==1 iff vertex v hasn't moved to the opposite group.
	 */
	unmoved = allocate_hash_set_unmoved(gSize);

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs, TRUE);
		maxi=0;
		maxScore = -DBL_MAX;

		k=getNextSetFlag(unmoved,gSize,k);
		while (k!=-1){
			Q_t = BEN_get_modularity_moved(Bg, s_ptr, Bs, k) - Q_0;
			if (Q_t > maxScore){
				maxScore = Q_t;
				maxi = k;
			}
			k=getNextSetFlag(unmoved, gSize, k);
		}
		/* line 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		modify_Bs(Bg, Bs, maxi , s_ptr[maxi]);

		/* line 14-18, 21 in Alg. 4 PsCode */
		tmpImprove += maxScore;
		if (tmpImprove > *maxImprove){
			*maxImprove = tmpImprove;
			*maxImpInd = i;
		}
		/* line 19 in Alg. 4 PsCode */
		resetFlag(unmoved, gSize, maxi);
	}
	free(unmoved);

	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: STEP 1 - move_maximal_score_vertex_mod_2\n");
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_mod_3(modMat *Bg, vector *s, vector Bs){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove=-DBL_MAX;
	/* clock_t start, end; */
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)

	
	#ifdef DEBUG_DIV_TWO
	srand(time(NULL));
	start = clock();
	printf("BEGIN: optimize_division_original_mod_3\n");
	#endif
	do {
		lastDeltaQ=deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex_mod_3(Bg, s, Bs, indices, &maxImprove, &maxImpInd);

		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j=indices+gSize-1; j>indices+maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		/* VERIFY(iter < gSize, INFINITE_LOOP_ERROR) */
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && ++iter < gSize);
	#ifdef DEBUG_DIV_TWO
	printf("SUCCESS: optimize_division_mod_3\n");
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif
	free(indices);
}

/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex_mod_Linked(modMat *Bg, vector *s, vector Bs, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	Linked_list_moved *unmoved = NULL;
	num i, maxi = 0, gSize = Bg->gSize;
	mNode *head, *prev = NULL, *delhead, *delprev = NULL;
	double Q_0, Q_t, maxScore = -DBL_MAX, tmpImprove=0;

	/* line 2 in Alg. 4 PsCode - a hash set of boolean values, s.t. 
	 * moved[v]==1 iff vertex v has moved to the opposite group.
	 */
	unmoved = create_unmoved(gSize);

	for (i = 0 ; i < gSize ; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		Q_0 = get_modularity_init(Bg, s_ptr, Bs, TRUE);
		maxi=0;
		maxScore = -DBL_MAX;
		head = unmoved->head;
		prev = NULL;
		delprev = NULL;
		delhead = head;
		while (head != NULL){
			Q_t = BEN_get_modularity_moved(Bg, s_ptr, Bs, head->ind) - Q_0;
			if (Q_t > maxScore){
				maxScore = Q_t;
				maxi = head->ind;
				delhead = head;
				delprev = prev;
			}
			prev = head;
			head = head->next;
		}
		
		/* line 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;
		modify_Bs(Bg, Bs, maxi , s_ptr[maxi]);

		/* line 14-18, 21 in Alg. 4 PsCode */
		tmpImprove += maxScore;
		if (tmpImprove > *maxImprove){
			*maxImprove = tmpImprove;
			*maxImpInd = i;
		}
		/* line 19 in Alg. 4 PsCode */
		delete_node(unmoved, delprev, delhead);
	}
	delete_unmoved(unmoved);
}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division_mod_Linked(modMat *Bg, vector *s, vector Bs){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove=-DBL_MAX;
	
	indices=(int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)


	do {
		lastDeltaQ=deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex_mod_Linked(Bg, s, Bs, indices, &maxImprove, &maxImpInd);

		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j=indices+gSize-1; j>indices+maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		/* VERIFY(iter < gSize, INFINITE_LOOP_ERROR) */
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)) && ++iter < gSize);

	#ifdef DEBUG_DIV_TWO
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
	clock_t start, end;

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
	
	srand(time(NULL));
	start = clock();
	printf("START: optimize_division_mod_3, on %d vertices\n", sizeG);

	optimize_division_mod_3(Bg, &s, u);

	end = clock();
	printf("END: optimize_division_mod_3, on %d vertices\n", sizeG);
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	Bg->free(Bg);

	map_s_to_groups(g, sizeG, s, g1, g2, sizeG1, sizeG2);

	free(u);
	free(s);
	
	return ret;
}