#include "Divide_Into_Two.h"

/** First part of Alg. 4 - find improvement of the partition defined by s */
void move_maximal_score_vertex(modMat *Bg, vector *s, int_vector indices, double *maxImprove, num *maxImpInd) {
	vector s_ptr = *s;
	long_num *unmoved=NULL;
	num i, maxi = 0, gSize = Bg->gSize;
	int k=0;
	double Q_t, maxScore = -DBL_MAX, tmpImprove=0, maxImprovement = -DBL_MAX;
	num maxImprovementInd = 0;
	

	/* initializing a flag set of boolean values, s.t. 
	 * unmoved[v]==1 iff vertex v hasn't moved to the opposite group.
	 */
	unmoved = allocate_flag_set(gSize);

	for (i = 0 ; i < gSize ; i++){
		maxi=0;
		maxScore = -DBL_MAX;
		k = get_next_set_flag(unmoved, gSize, 0, TRUE);
		while (k != -1){
			Q_t = Bg->get_modularity(Bg, s_ptr, k);
			if (Q_t > maxScore){ /* saving the maximum score achieved and the index j' caused it */
				maxScore = Q_t;
				maxi = k;
			}
			k = get_next_set_flag(unmoved, gSize, k, FALSE);
		}
		s_ptr[maxi] *= -1;
		*(indices++)  = maxi;

		tmpImprove += maxScore; /* updating improve[i] with score[j'] */
		if (tmpImprove > maxImprovement){ /* saving the maximum improvement */
			maxImprovement = tmpImprove;
			maxImprovementInd = i;
		}
		reset_flag(unmoved, gSize, maxi); /* Unmoved = Unmoved \ {j'} */
	}
	free(unmoved);
	*maxImprove = maxImprovement;
	*maxImpInd = maxImprovementInd;
	
}


/** Optimizes a 2-division by moving vertices in one group to the opposite group and attempting to ascend modularity Q of Bg.
 * 	The function will modify vector s, if an optimized division is found.
 * 	Based on Alg. 4 pseudo-code.
 * 	@param Bg - The modularity matrix of a network.
 * 	@param s - a {-1,1} vector encoding a 2-division.
 **/
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
		/* find improvement of the partition defined by s */
		move_maximal_score_vertex(Bg, s, indices, &maxImprove, &maxImpInd);
		/* find the maximum improvement of s and update s accordingly */
		for (j= indices + gSize - 1; j > indices + maxImpInd; j--)
			s_ptr[*j] *= -1;
		deltaQ = (maxImpInd == gSize-1) ? 0 : maxImprove;

		VERIFY(iter++ < maximize_loop_limit, INFINITE_LOOP_ERROR)
	} while (IS_POSITIVE(deltaQ) && IS_POSITIVE(fabs(deltaQ-lastDeltaQ)));
	
	free(indices);
}


/** Maps a real vector to a pre-allocated vector of {-1,1} elements. */
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

/** Maps a {-1,1} vector s s.t. dim(s)==dim(B) to a partition g1, g2 of g,
 *  where g is a subset of 0,1,...,n, s.t. n=|V| of the original network.
 **/  
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

/* Algorithm 2 */
DIV_RESULT divide_into_two(modMat *B, Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2){
	scalar beta;
	vector u, s;
	modMat *Bg;
	DIV_RESULT ret;

	/* creating a sub matrix B_hat[g] */
	Bg = create_Sub_Matrix(B, g, sizeG);

	/* computing leading eigenpair u1 and B1 of the modularity matrix B_hat[g] */
	beta = Leading_eigenpair(Bg, &u);

	if (!IS_POSITIVE(beta)){
		Bg->free(Bg);
		free(u);
		*sizeG1 = sizeG;
		*sizeG2 = 0;
		return GROUP_INDIVISIBLE;
	}
	/* computing s according to u1 */
	eigen_to_s(Bg, u, &s);

	if (!IS_POSITIVE(Bg->get_modularity(Bg, s, MODULARITY_INIT)))
		ret=GROUP_INDIVISIBLE;
	else
		ret=GROUP_DIVIDED;
	optimize_division(Bg, &s);

	/* B_hat[g] is no longer in use */
	Bg->free(Bg);

	/* updating g1 and g2 for algorithm 3 */
	map_s_to_groups(g, sizeG, s, g1, g2, sizeG1, sizeG2);

	free(u);
	free(s);

	return ret;
}