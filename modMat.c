#include "modMat.h"

/* multiply v*u */
double dot_prod(vector v, vector u, num d){
	num i;
	double acc=0;
	for (i=0; i<d; i++, v++, u++){
		acc += (*v)*(*u);
	}
	return acc;
}


/* helping function, populating Bg->K, Bg->currM corresponding to the subgroup g*/
void compute_K_and_currM(modMat *Bsrc, modMat *Bg, Subgroup g, num gSize){
	int_vector g_i, K_i=Bsrc->K,  Kg_i=Bg->K;
    num cnt = 0;
    for (g_i = g, Kg_i = Bg->K; g_i < gSize + g ; g_i++, Kg_i++){
        *Kg_i = K_i[*g_i];
        cnt += *Kg_i;
    }
	Bg->currM=cnt;
}

/* Constructor, creates a new sub matrix B_hat[g], corresponds to B and g */
modMat* create_Sub_Matrix(modMat *B, Subgroup g, num sizeG){
    modMat *Bg=NULL;
	Bg = allocate_mod_mat(sizeG, 0, TRUE);
	VERIFY(Bg != NULL, MEM_ALLOC_ERROR)
	Bg->A = B->A->create_sub_mat(B->A, g, sizeG); /* create a sub spmat */
	VERIFY(Bg->A != NULL, MEM_ALLOC_ERROR)
	compute_K_and_currM(B, Bg, g, sizeG);
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
	Bg->one_norm = B->one_norm; /* Use the greatest 1-norm for all matrices of any size less than B->gSize */
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) minus (f_i - ||C||) * I (while ||C|| is the one norm of Bg, in case a shift is needed),
 * all multiplications at one loop
 * */
void mult_K_F_and_C(const modMat *Bg, const double *v, double *res, boolean shift){
    int_vector K = Bg->K, spmatSize = Bg->A->spmatSize;
	num currM = Bg->currM ,sizeG = Bg->gSize, *ki;
	double origM = (double)Bg->M; 
    double dot = 0;
	double fi, shiftNorm, KFC;
	
	VERIFY(IS_POSITIVE(origM), DIVISION_BY_ZERO)
    VERIFY(res != NULL,NULL_POINTER_ERROR)

    for (ki = K ; ki < sizeG + K ; ki++, v++){
        dot += (*ki) * (*v);
    }

	v -= sizeG;
	shiftNorm = shift ? Bg->one_norm : 0;

    for (ki = K ; ki < sizeG + K ; ki++, v++, res++, spmatSize++){
		KFC = 0; /* KFC stores the multiplication of the degrees matrix (k_i * k_j / M) minus (f_i - ||C||) * I --- i.e. (K + f_i - ||C||)*v*/
		if (*ki != 0){ /* an integer value non-zero check */
			KFC = (*ki) * dot * (1 / origM);
			fi = (*spmatSize) - ((currM * (*ki)) / origM);
			KFC += (fi * (*v));
		}
		KFC -= (shiftNorm * (*v));
		*res -= KFC; /* res currently stores Av, so the overall result will be = (Av - (K + f_i - ||C||)*v) */
    }
	
}


/* Implements multiplication of B_hat with a vector by
 * using several mult. functions and adding results together.
 */
void mult_B_hat(const struct _modmat *Bg, const double *v, double *result, boolean shift){	
	Bg->A->mult(Bg->A, v, result);
	mult_K_F_and_C(Bg, v, result, shift);
}

double sum_of_abs(double *row, num n){
	double *p, res=0;
	for (p = row ; p < row + n ; p++)
		res += fabs(*p);
	return res;
}

void set_1_norm(modMat *B){
	/*
	 * By symmetry, row i == column i.
	 * Thus, any operation mapped (entry-wise) on a row of B is equivalent to the same operation on a column of B.
	 */
	num i, gSize = B->gSize;
	double tmp = 0, max = 0;
	vector B_i, e_i, e_i_ptr;
	
	B_i = (vector)malloc(gSize * sizeof(double));
	VERIFY(B_i!=NULL, MEM_ALLOC_ERROR)

	e_i = (vector)calloc(gSize,sizeof(double));
	VERIFY(e_i!=NULL,MEM_ALLOC_ERROR)
	e_i_ptr = e_i;

	*e_i_ptr = 1;
	for (i = 0 ; i < gSize-1 ; i++){
		B->mult(B, e_i, B_i, NO_SHIFT);
		tmp = sum_of_abs(B_i, gSize);
		if (tmp > max)
			max = tmp;
		*e_i_ptr++ = 0;
		*e_i_ptr = 1;
	}

	/* Last row of B */
	B->mult(B, e_i, B_i, NO_SHIFT);
	tmp = sum_of_abs(B_i, gSize);
	if (tmp > max)
		max = tmp;

	free(B_i);
	free(e_i);
	B->one_norm = max;
}

/** Compute Modularity of B[g]_hat: 0.5 * s^T * B[g]_hat * s.
** the functions can be called from the Maximization algorithm, where only one entry of s is being changed, so the spmat module can calculate it smartly and quickly.
** @param move_vertex - if MODULARITY_INIT, compute modularity w.r.t to S. Otherwise, move vertex s[move_vertex] temporarily and compute.
*/
double get_B_modularity(struct _modmat *B, vector s, int move_vertex){
	double dQ;
	vector Bs;
	if (move_vertex == MODULARITY_INIT){
		Bs = (vector)malloc(B->gSize*sizeof(double));
		VERIFY(Bs!=NULL, MEM_ALLOC_ERROR)

		B->mult(B, s, Bs, NO_SHIFT);
		dQ = dot_prod(s, Bs, B->gSize);
		
		free(Bs);
	}
	else
		dQ = B->A->get_modularity_score(B->A, s, move_vertex, B->K, B->M);
	return dQ;
}

/* free B's allocated fields and B itselft */
void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B);
}

/* allocation function, can allocate both sub modMat and full size modMat */
modMat* allocate_mod_mat(num n, num nnz, boolean isSub){
	modMat *rep = (modMat*)malloc(sizeof(modMat));
	VERIFY(rep!=NULL, MEM_ALLOC_ERROR)

	rep->gSize = n;
	if (!(isSub)){
		if (USE_SPMAT_LINKED)
			rep->A = spmat_allocate_list(n);
		else
			rep->A = spmat_allocate_array(n, nnz);
		VERIFY(rep->A != NULL,MEM_ALLOC_ERROR)
	}
	rep->K = (int_vector)malloc(n * sizeof(num));
	VERIFY(rep->K != NULL,MEM_ALLOC_ERROR)

	rep->free=free_mod_mat;
	rep->mult=mult_B_hat;
	rep->get_modularity=get_B_modularity;

	return rep;
}
