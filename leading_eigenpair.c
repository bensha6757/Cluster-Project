#include "Leading_eigenpair.h"

double l2_norm(vector v, num d){
	return sqrt(dot_prod(v,v,d));
}

void set_rand_vector(vector v, num n){
	vector p;
	for (p=v; p<v+n; p++)
		*p = (double)rand();
}

double l1_dist(vector a, vector b, num len){
	double maxDiff=0, diff;
	vector p;
	for (p=a; p < a+len; p++, b++){
		diff = fabs(*p - *b);
		if (diff > maxDiff)
			maxDiff = diff;
	}
	return maxDiff;
}

/*
 * Apply a single iteration to a matrix B and a vector v according to the power method.
 * Store result in result vector.
 */
void power_iteration(modMat *Bg, vector v, vector result){
	vector p;
	double nrm;
	num n=Bg->gSize;
	Bg->mult(Bg, v, result, SHIFT);
	nrm = l2_norm(result, n);
	VERIFY(IS_POSITIVE(nrm), DIVISION_BY_ZERO)
	for (p = result; p < result + n; p++)
		*p /= nrm;
}


/**	Approximate dominant eigen value of matrix B, 
 * 	using vectors computed in last power iteration, according to equation:
 *  beta_1 = (Ab_k * b_k) / ||b_k||^2.
 **/
double approx_dom_eigen_val(modMat *B, vector Ab_k, vector b_k){
	num gSize=B->gSize;
	double d;
	B->mult(B, b_k, Ab_k, SHIFT);
	d = dot_prod(b_k, b_k ,gSize);
	VERIFY(IS_POSITIVE(d),DIVISION_BY_ZERO)
	return dot_prod(b_k, Ab_k ,gSize)/d;
}

scalar Leading_eigenpair(modMat *Bg, vector *leadEigenVec){
	num iter=0, gSize = Bg->gSize;
	num eigenpair_loops_bound = ITER_LIMIT(gSize);
	vector bprev, bnext, tmp;
	scalar leadEigenValue;
	scalar last_l1_dist = 0;
	
	bprev = (vector)malloc(gSize * sizeof(scalar));
	VERIFY(bprev!=NULL,MEM_ALLOC_ERROR)
	bnext = (vector)malloc(gSize * sizeof(scalar));
	VERIFY(bnext!=NULL,MEM_ALLOC_ERROR)

	set_rand_vector(bnext, gSize);

	do {

		tmp = bnext;
		bnext = bprev;
		bprev = tmp;
		power_iteration(Bg, bprev, bnext);

		last_l1_dist = l1_dist(bprev, bnext, gSize);
		VERIFY(iter++ < eigenpair_loops_bound, INFINITE_LOOP_ERROR)

	} while ( IS_POSITIVE(last_l1_dist) && iter < eigenpair_loops_bound  );

	*leadEigenVec = bnext;

	/* The leading eigenvalue is a 1-norm shifted dominant eigenvalue*/
	leadEigenValue = approx_dom_eigen_val(Bg,bprev,bnext) - (Bg->one_norm);
	
	free(bprev);
	return leadEigenValue;
}