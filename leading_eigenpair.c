#include "leading_eigenpair.h"

double dot_prod(vector v, vector u, num d){
	num i;
	double acc=0;
	for (i=0; i<d; i++, v++, u++){
		acc += (*v)*(*u);
	}
	return acc;
}

double l2_norm(vector v, num d){
	return sqrt(dot_prod(v,v,d));
}

void set_rand_vector(vector v, num n){
	vector p;
	for (p=v; p<v+n; p++)
		*p = (double)rand();

	#ifdef DEBUG_EIGEN
	printf("SUCCESS: set_rand_vector\n");
	#endif
}

/*	Approximate dominant eigen value of matrix B, using vectors computed in last power iteration */
double approx_dom_eigen_val(modMat *B, vector bprev, vector bnext){
	num gSize=B->gSize;
	double d;
	B->mult(B,bnext,bprev,SHIFT);
	d = l2_norm(bnext,gSize);
	/*return dot_prod(bnext,bprev,gSize) / dot_prod(bprev,bprev,gSize);*/
	return dot_prod(bnext,bprev,gSize)/(d*d);
}

boolean is_within(vector a, vector b, num d){
	num i;
	for (i=0; i<d; i++, a++, b++){
		if (IS_POSITIVE(fabs(*a - *b)))
			return 0;
	}
	return 1;
}

/*
 * Apply a single iteration to a matrix B and a vector v according to the power method.
 * Store result in result vector.
 */
void power_iteration(modMat *Bg, vector v, vector result){
	vector p;
	double nrm;
	num n=Bg->gSize;
	#ifdef DEBUG_EIGEN
	printf("BEGIN: power_iteration\n");
	#endif
	Bg->mult(Bg, v, result, SHIFT);
	nrm = l2_norm(result, n);
	for (p = result; p < result + n; p++)
		*p /= nrm;
	#ifdef DEBUG_EIGEN
	printf("SUCCESS: power_iteration\n");
	#endif
}

/*
 * Compute leading eigen pair of modularity Matrix B_hat[g].
 */
void leading_eigenpair(modMat *Bg, vector *leadEigenVec, scalar *leadEigenVal){
	num iter=0, gSize = Bg->gSize;
	num loops_limit = (num)pow(gSize,2); /*There are about O(N^2) Power iterations for a matrix of size N */
	vector bprev, bnext;
	/*vector p, q;*/
	
	#ifdef DEBUG_EIGEN
	printf("BEGIN: leading_eigenpair\n");
	#endif
	bprev = (vector)malloc(gSize * sizeof(double));
	VERIFY(bprev!=NULL,MEM_ALLOC_ERROR)
	bnext = (vector)malloc(gSize * sizeof(double));
	VERIFY(bnext!=NULL,MEM_ALLOC_ERROR)

	srand(time(NULL));
	set_rand_vector(bnext, gSize);

	do {
		memcpy(bprev, bnext, gSize * sizeof(double));
		power_iteration(Bg, bprev, bnext);
		iter++;
		/*VERIFY(iter < 1000 * gSize,INFINITE_LOOP_ERROR)*/
	} while (!is_within(bprev, bnext, gSize) && iter < loops_limit);
	#ifdef DEBUG_EIGEN
		printf("# of power iterations: %d\n", (int)iter);
	#endif
	*leadEigenVec = (vector)malloc(gSize * sizeof(double));
	VERIFY(*leadEigenVec!=NULL, MEM_ALLOC_ERROR)
	memcpy(*leadEigenVec,bnext,gSize * sizeof(double));

	/* The leading eigenvalue is a 1-norm shifted dominant eigenvalue*/
	*leadEigenVal = approx_dom_eigen_val(Bg,bprev,bnext) - (Bg->one_norm);
	free(bnext);
	free(bprev);
	#ifdef DEBUG_EIGEN
	printf("SUCCESS: leading_eigenpair, %f\n", *leadEigenVal);
	#endif
}