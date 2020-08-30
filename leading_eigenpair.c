#include "leading_eigenpair.h"
#define DEBUG

double dot_prod(vector v, vector u, num d){
	unsigned int i;
	double acc=0;
	for (i=0; i<d; i++, v++, u++){
		acc += (*v)*(*u);
	}
	return acc;
}

double norm(vector v, num d){
	return sqrt(dot_prod(v,v,d));
}

void set_rand_vector(vector v, num n){
	unsigned int i;
	srand(time(NULL));
	for (i=0; i<n; i++)
		*v++ = (double)rand();
}

/*	Approximate dominant eigen value of matrix B, using vectors computed in last power iteration */
double approx_dom_eigen_val(modMat *B, vector bprev, vector bnext){
	return dot_prod(bnext,bprev,B->gSize) / dot_prod(bprev,bprev,B->gSize);
}

unsigned int is_within(vector a, vector b, num d){
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
	mult_B_hat_g(Bg, v, result, SHIFT);
	nrm = norm(result, Bg->gSize);
	for (p = result; p < result + Bg->gSize; p++)
		*p /= nrm;
}

/*
 * Compute leading eigen pair of modularity Matrix B_hat[g].
 */
void leading_eigenpair(modMat *Bg, vector *leadEigenVec, scalar *leadEigenVal){
	num iter=0, gSize = Bg->gSize;
	num loops_limit = gSize*gSize*gSize;
	vector bprev, bnext, tmp;
	
	#ifdef DEBUG
	printf("BEGIN: leading_eigenpair\n");
	#endif

	*leadEigenVec = (vector)malloc(gSize * sizeof(double));
	VERIFY(*leadEigenVec!=NULL, MEM_ALLOC_ERROR)
	bprev = (vector)malloc(gSize * sizeof(double));
	VERIFY(bprev!=NULL,MEM_ALLOC_ERROR)
	bnext = (vector)malloc(gSize * sizeof(double));
	VERIFY(bnext!=NULL,MEM_ALLOC_ERROR)
	set_rand_vector(bprev, gSize);
	power_iteration(Bg, bprev, bnext);
	iter++;
	while (!is_within(bprev, bnext, gSize) && iter < loops_limit){
		memcpy(bprev, bnext, gSize * sizeof(double));
		power_iteration(Bg, bprev, bnext);
		iter++;
		/*VERIFY(iter < 1000 * gSize,INFINITE_LOOP_ERROR)*/
	}
	#ifdef DEBUG
		printf("# of power iterations: %d\n", (int)iter);
	#endif
	tmp = *leadEigenVec;
	memcpy(tmp,bnext,gSize * sizeof(double));
	/* The leading eigenvalue is a 1-norm shifted dominant eigenvalue*/
	*leadEigenVal = approx_dom_eigen_val(Bg,bprev,bnext) - (Bg->one_norm);
	free(bnext);
	free(bprev);
	#ifdef DEBUG
	printf("SUCCESS: leading_eigenpair\n");
	#endif
}