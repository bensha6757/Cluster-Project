#include "modMat.h"

void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->spmatSize);
	free(B);
}


void get_K_row(const modMat *Bg, num i, double *row){
	int_vector p, K = Bg->K;
	num k_i, gSize = Bg->gSize;

	k_i = K[i];
	for (p = K; p < K + gSize; p++, row++){
		*row = (k_i) * ((*p) / (double)Bg->M);
	}
}

/* Copy row i of B_hat matrix of size B->gSize to row vector */

void get_B_hat_row(const struct _modmat *Bg, num i, double *row){
	double *A_i, *K_i, *r;
	double f_i = 0;
	num gSize = Bg->gSize;
	#ifdef DEBUG_MODMAT
	printf("BEGIN: get_B_hat_row=%d\n", i);
	#endif
	A_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(A_i!=NULL, MEM_ALLOC_ERROR)
	K_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(K_i!=NULL, MEM_ALLOC_ERROR)
	Bg->A->get_row(Bg->A, i, A_i);
	get_K_row(Bg, i, K_i);

	for (r=row; r < row + gSize; r++, A_i++, K_i++){
		*r = *A_i - *K_i;
		f_i += *r;
	}
	row[i] -= f_i;
	free(A_i - gSize);
	free(K_i - gSize);
}


void get_B_hat_row_generic(const struct _modmat *Bg, num i, double *row){
	vector e_i;
	get_basis_unit_vec(&e_i, i, Bg->gSize);
	Bg->mult(Bg, e_i, row, NO_SHIFT);
	free(e_i);
}

double sum_of_abs(double *row, num n){
	double *p, res=0;
	for (p = row ; p < row + n ; p++)
		res += fabs(*p);
	return res;
}

/* helping function, populating Bg->K, Bg->currM according to the subgroup g*/
void compute_K_and_currM(modMat *Bg, int_vector K, Subgroup g, num gSize){
	int_vector g_i, Kg_i;
    num Mg = 0;

    for (g_i = g, Kg_i = Bg->K; Kg_i < gSize + Bg->K ; g_i++, Kg_i++){
        *Kg_i = K[*g_i];
        Mg += *Kg_i;
    }
	
    Bg->currM = Mg;
}



/* Constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */

modMat *create_Sub_Matrix(modMat *B, Subgroup g, num sizeG){
    modMat *Bg;

	Bg = allocate_mod_mat(sizeG, sizeG*sizeG);
    VERIFY(Bg!=NULL,MEM_ALLOC_ERROR)
	if (USE_SPMAT_LINKED)
		Bg->A = create_sub_sparse_matrix_linked(B->A, g, sizeG, Bg->spmatSize);
	else
		Bg->A = create_sub_sparse_matrix_array(B->A, g, sizeG, Bg->spmatSize);
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
    compute_K_and_currM(Bg, B->K, g, sizeG);
	/*set_1_norm(Bg);*/
	Bg->one_norm = B->one_norm; /* Use the greatest 1-norm for all matrices of size < B->gSize */
	
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) 
 * */
void mult_K(const modMat *Bg, const double *v, double *res){
    int_vector K = Bg->K;
    num origM = Bg->M, sizeG = Bg->gSize, *ki;
    double dot = 0;
	
    VERIFY(res != NULL,NULL_POINTER_ERROR)
    for (ki = K ; ki < sizeG + K ; ki++, v++){
        dot += (*ki) * (*v);
    }
    for (ki = K ; ki < sizeG + K ; ki++, res++){
      *res = (*ki) * dot * (1 / (double)origM);
    }
	
}


/** part of the modularity matrix multiplication for power iteration, 
 * 	multiplying the 2 matrices (f_i - ||C||) * I 
 * */
void mult_F_and_C(const modMat *Bg, const double *v, double *res, boolean shift){
    int_vector K = Bg->K;
    int_vector spmatSize = Bg->spmatSize;
    num M = Bg->currM, origM = Bg->M ,sizeG = Bg->gSize, *ki;
    double fi, shiftNorm;
	
	shiftNorm = shift ? Bg->one_norm : 0;
    for (ki = K ; ki < sizeG + K ; ki++, v++, res++, spmatSize++){
		fi = (*spmatSize) - (((*ki) * M) / (double)origM);
		*res = (fi - shiftNorm)  * (*v);
    }
	
}


/* Implements multiplication of B_hat with a vector by
 * using several mult. functions and adding results together.
 */
void mult_B_hat(const struct _modmat *Bg, const double *v, double *result, boolean shift){
	vector Kv, Fv_minus_Cv, p;
	num gSize = Bg->gSize;
	
	Bg->A->mult(Bg->A,v,result);
	
	Kv=(vector)malloc(sizeof(double)*gSize);
	VERIFY(Kv!=NULL,MEM_ALLOC_ERROR)
	mult_K(Bg, v, Kv);

	Fv_minus_Cv=(vector)malloc(sizeof(double)*gSize);
	VERIFY(Fv_minus_Cv!=NULL,MEM_ALLOC_ERROR)
	mult_F_and_C(Bg, v, Fv_minus_Cv, shift);

	for (p=result; p < result + gSize; p++, Kv++, Fv_minus_Cv++){
		*p -= *Kv;
		*p -= *Fv_minus_Cv;
	}
	
	free(Kv-gSize);
	free(Fv_minus_Cv-gSize);	
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
		B->get_row(B,i,B_i);
		tmp = sum_of_abs(B_i, gSize);
		if (tmp > max)
			max = tmp;
	}
	free(B_i);
	B->one_norm = max;
}

/*allocate new ModMat of size n*/
modMat* allocate_mod_mat(num n, num nnz){
	modMat *rep = (modMat*)malloc(sizeof(modMat));
	VERIFY(rep!=NULL,MEM_ALLOC_ERROR)
	rep->gSize = n;
	if (USE_SPMAT_LINKED)
		rep->A = spmat_allocate_list(n);
	else
		rep->A = spmat_allocate_array(n,nnz);
	VERIFY(rep->A != NULL,MEM_ALLOC_ERROR)

	rep->K = (int_vector)malloc(n * sizeof(num));
	VERIFY(rep->K != NULL,MEM_ALLOC_ERROR)
	
	rep->spmatSize=(int_vector)malloc(n * sizeof(num));
	VERIFY(rep->spmatSize != NULL,MEM_ALLOC_ERROR)

	rep->free=free_mod_mat;
	rep->mult=mult_B_hat;
	rep->get_row=get_B_hat_row;

	return rep;
}
