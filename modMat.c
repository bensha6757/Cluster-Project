#include "modMat.h"

double dot_prod(vector v, vector u, num d){
	num i;
	double acc=0;
	for (i=0; i<d; i++, v++, u++){
		acc += (*v)*(*u);
	}
	return acc;
}


/* helping function, populating Bg->K, Bg->currM according to the subgroup g*/
void compute_K_and_currM(modMat *Bsrc, modMat *Bg, Subgroup g, num gSize){
	int_vector g_i, K_i=Bsrc->K,  Kg_i=Bg->K;
    num cnt = 0;
    for (g_i = g, Kg_i = Bg->K; g_i < gSize + g ; g_i++, Kg_i++){
        *Kg_i = K_i[*g_i];
        cnt += *Kg_i;
    }
	Bg->currM=cnt;
}

modMat* create_Sub_Matrix(modMat *B, Subgroup g, num sizeG){
    modMat *Bg=NULL;
	spmat *Ag;
	Bg = allocate_mod_mat(sizeG, 0, TRUE);
	VERIFY(Bg != NULL, MEM_ALLOC_ERROR)
	Ag = B->A->create_sub_mat(B->A, g, sizeG);
	VERIFY(Ag != NULL, MEM_ALLOC_ERROR)
	compute_K_and_currM(B, Bg, g, sizeG);
	Bg->A = Ag;
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
	Bg->one_norm = B->one_norm; /* Use the greatest 1-norm for all matrices of any size less than B->gSize */
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) minus (f_i - ||C||) * I
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
		KFC = 0;
		if (*ki != 0){ /* an integer value non-zero check */
			KFC = (*ki) * dot * (1 / origM);
			fi = (*spmatSize) - ((currM * (*ki)) / origM);
			KFC += (fi * (*v));
		}
		KFC -= (shiftNorm * (*v));
		*res -= KFC;
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


double get_B_modularity_score(modMat *Bg, vector s, int moved_v){
	double dQ;
	vector A_j, d_j;
	int_vector K = Bg->K, K_j;
	num gSize=Bg->gSize, k_i = K[moved_v];
	int d_i;
	double k_i_M = k_i / (double) Bg->M, sum = 0;

	A_j = (vector)malloc(gSize*sizeof(double));
	VERIFY(A_j!=NULL, MEM_ALLOC_ERROR)
	Bg->A->get_row(Bg->A, moved_v, A_j);
	
	s[moved_v] *= -1;
	d_i = s[moved_v];

	for(d_j = s, K_j = K; d_j < gSize + s ; d_j++, K_j++, A_j++){
		sum += ((*A_j - (*K_j * k_i_M)) * (*d_j));
	}

	dQ = (4 * d_i * sum) + (4 * k_i_M * k_i);
	/* Restore s to initial state */
	s[moved_v] *= -1;

	A_j -= gSize;
	free(A_j);
	
	return dQ;
}

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


void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B);
}


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
