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

/* Constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */

modMat* create_Sub_Matrix(modMat *B, Subgroup g, num sizeG){
    modMat *Bg=NULL;
	spmat *Ag;
	num nnz=0;
	Ag = B->A->create_sub_mat(B->A, g, sizeG);
	VERIFY(Ag!=NULL, MEM_ALLOC_ERROR)
	if (!USE_SPMAT_LINKED)
		nnz=sum_array(Ag->spmatSize, sizeG);
	Bg = allocate_mod_mat(sizeG, nnz, TRUE);
	VERIFY(Bg!=NULL, MEM_ALLOC_ERROR)
	compute_K_and_currM(B, Bg, g, sizeG);
	Bg->A = Ag;
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
	Bg->one_norm = B->one_norm; /* Use the greatest 1-norm for all matrices of size < B->gSize */
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) minus (f_i - ||C||) * I
 * */
void mult_K_F_and_C(const modMat *Bg, const double *v, double *res, boolean shift){
    int_vector K = Bg->K, spmatSize = Bg->A->spmatSize;
	num currM = Bg->currM, origM = Bg->M ,sizeG = Bg->gSize, *ki;
    double dot = 0;
	double fi, shiftNorm, KFC;
	
    VERIFY(res != NULL,NULL_POINTER_ERROR)
    for (ki = K ; ki < sizeG + K ; ki++, v++){
        dot += (*ki) * (*v);
    }
	v -= sizeG;
	shiftNorm = shift ? Bg->one_norm : 0;

    for (ki = K ; ki < sizeG + K ; ki++, v++, res++, spmatSize++){
		KFC = 0;
		if (*ki != 0){
			KFC = (*ki) * dot * (1 / (double)origM);
			fi = (*spmatSize) - ((currM * (*ki)) / (double)origM);
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

double get_K_and_F_row(const struct _modmat *Bg, num i, double *row){
	double *r, f_i = 0;
	int_vector k_j;
	num gSize = Bg->gSize, k_i = Bg->K[i];
	
	for (r = row, k_j = Bg->K; r < row + gSize; r++, k_j++){
		*r -= ((k_i) * ((*k_j) / (double)Bg->M));
		f_i += *r;
	}

	return f_i;
}

/* Copy row i of B_hat matrix of size B->gSize to row vector */
void get_B_hat_row(const struct _modmat *Bg, num i, double *row){
	double f_i;
	Bg->A->get_row(Bg->A, i, row);
	f_i = get_K_and_F_row(Bg, i, row);
	row[i] -= f_i;
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
	vector B_i;
	
	B_i = (vector)malloc(gSize * sizeof(double));
	VERIFY(B_i!=NULL,MEM_ALLOC_ERROR)
	for (i = 0 ; i < gSize ; i++){
		B->get_row(B, i, B_i);
		tmp = sum_of_abs(B_i, gSize);
		if (tmp > max)
			max = tmp;
	}
	free(B_i);
	B->one_norm = max;
}



void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B);
}

double get_B_modularity(struct _modmat *B, vector s, vector Bs){
	double Q;
	
	B->mult(B, s, Bs, NO_SHIFT);
	Q = dot_prod(s, Bs, B->gSize);
	
	return Q;
}



/*allocate new ModMat of size n*/
modMat* allocate_mod_mat(num n, num nnz, boolean isSub){
	modMat *rep = (modMat*)malloc(sizeof(modMat));
	VERIFY(rep!=NULL,MEM_ALLOC_ERROR)

	rep->gSize = n;

	if (!(isSub)){
		if (USE_SPMAT_LINKED)
			rep->A = spmat_allocate_list(n);
		else
			rep->A = spmat_allocate_array(n,nnz);
		VERIFY(rep->A != NULL,MEM_ALLOC_ERROR)
	}

	rep->K = (int_vector)malloc(n * sizeof(num));
	VERIFY(rep->K != NULL,MEM_ALLOC_ERROR)

	rep->free=free_mod_mat;
	rep->mult=mult_B_hat;
	rep->get_row=get_B_hat_row;
	rep->get_modularity=get_B_modularity;

	return rep;
}
