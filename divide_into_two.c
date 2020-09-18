#include "Divide_Into_Two.h"

/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */
void move_maximal_score_vertex(modMat *Bg, vector *s, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	long_num *unmoved=NULL;
	num i, maxi = 0, gSize = Bg->gSize;
	int k=0;
	double Q_t, maxScore = -DBL_MAX, tmpImprove=0, maxImprovement = -DBL_MAX;
	num maxImprovementInd = 0;
	

	/* line 2 in Alg. 4 PsCode - a flag set of boolean values, s.t. 
	 * unmoved[v]==1 iff vertex v hasn't moved to the opposite group.
	 */
	unmoved = allocate_flag_set(gSize);

	for (i = 0 ; i < gSize ; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		maxi=0;
		maxScore = -DBL_MAX;
		k = get_next_set_flag(unmoved, gSize, 0, TRUE);
		while (k != -1){
			Q_t = Bg->get_modularity(Bg, s_ptr, k);
			if (Q_t > maxScore){
				maxScore = Q_t;
				maxi = k;
			}
			k = get_next_set_flag(unmoved, gSize, k, FALSE);
		}
		/* line 12-13 in Alg. 4 PsCode */
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;

		/* line 14-18, 21 in Alg. 4 PsCode */
		tmpImprove += maxScore;
		if (tmpImprove > maxImprovement){
			maxImprovement = tmpImprove;
			maxImprovementInd = i;
		}
		/* line 19 in Alg. 4 PsCode */
		reset_flag(unmoved, gSize, maxi);
	}
	free(unmoved);
	*maxImprove = maxImprovement;
	*maxImpInd = maxImprovementInd;
	
}


/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */
void optimize_division(modMat *Bg, vector *s){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	num maximize_loop_limit = gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove=-DBL_MAX;

	indices = (int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)

	do {
		lastDeltaQ = deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex(Bg, s, indices, &maxImprove, &maxImpInd);

		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j= indices + gSize - 1; j > indices + maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		VERIFY(iter++ < maximize_loop_limit, INFINITE_LOOP_ERROR)
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)));
	
	free(indices);
}


/** Based on lines 2-20 in Alg. 4 pseudo-code .
 */


void move_maximal_score_vertex_mod_Linked(modMat *Bg, vector *s, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	Linked_list_moved *unmoved = NULL;
	num i, maxi = 0, gSize = Bg->gSize;
	mNode *head, *prev = NULL, *delhead, *delprev = NULL;
	double Q_t, maxScore = -DBL_MAX, tmpImprove=0, maxImprovement = -DBL_MAX;
	num maxImprovementInd = 0;

	unmoved = create_unmoved(gSize);
	head = unmoved->head;

	for (i=0; i<gSize; i++){
		/* lines 3-11 in Alg. 4 PsCode */
		maxi=0;
		maxScore = -DBL_MAX;
	
		while (head != NULL){
			Q_t = Bg->get_modularity(Bg, s_ptr, head->ind);
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

		/* line 14-18, 21 in Alg. 4 PsCode */
		tmpImprove += maxScore;
		if (tmpImprove > maxImprovement){
			maxImprovement = tmpImprove;
			maxImprovementInd = i;
		}
		/* line 19 in Alg. 4 PsCode */
		delete_node(unmoved, delprev, delhead);

		head = unmoved->head;
		prev = NULL;
	}
	delete_unmoved(unmoved);
	*maxImprove = maxImprovement;
	*maxImpInd = maxImprovementInd;

}



/** Optimize a division encoded by {-1,1} vector s by moving a vertex to other group and ascending modularity Q.
 * Based on lines 2-20 in Alg. 4 pseudo-code .
 * */

void optimize_division_mod_Linked(modMat *Bg, vector *s){
	vector s_ptr = *s;
	int_vector indices, j;
	num maxImpInd=0, iter=0, gSize = Bg->gSize;
	num maximize_loop_limit = gSize;
	double deltaQ=0, lastDeltaQ=0, maxImprove=-DBL_MAX;

	indices = (int_vector)malloc(gSize*sizeof(num));
	VERIFY(indices!=NULL, MEM_ALLOC_ERROR)

	do {
		lastDeltaQ = deltaQ;
		/* lines 2-21 in Alg. 4 pseudo-code */
		move_maximal_score_vertex_mod_Linked(Bg, s, indices, &maxImprove, &maxImpInd);
		/* lines 22-25 in Alg. 4 pseudo-code */
		for (j= indices + gSize - 1; j > indices + maxImpInd; j--)
			s_ptr[*j] *= -1;
		/* lines 26-30 in Alg. 4 pseudo-code */
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		VERIFY(iter++ < maximize_loop_limit, INFINITE_LOOP_ERROR)
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)));

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

	beta = Leading_eigenpair(Bg, &u);

	if (!IS_POSITIVE(beta)){
		/*printf("no optimize for size %d\n", sizeG);*/
		Bg->free(Bg);
		free(u);
		*sizeG1 = sizeG;
		*sizeG2 = 0;
		return GROUP_INDIVISIBLE;
	}
	
	eigen_to_s(Bg, u, &s);
	if (!IS_POSITIVE(Bg->get_modularity(Bg, s, MODULARITY_INIT)))
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;
	optimize_division(Bg, &s);

	Bg->free(Bg);

	map_s_to_groups(g, sizeG, s, g1, g2, sizeG1, sizeG2);

	free(u);
	free(s);

	return ret;
}