#include "modMat.h"

void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->g);
	free(B->spmatSize);
	free(B);
}


void get_K_row(const modMat *B, num i, double *row){
	int_vector p;
	num k_i;
	#ifdef DEBUG
	printf("BEGIN: get_K_row=%d\n", i);
	#endif
	k_i=*(B->K+i);
	for (p = B->K; p < B->K + B->gSize; p++)
		*row++=(k_i)*((double)((*p)/B->M));
	#ifdef DEBUG
	printf("SUCCESS: get_K_row=%d\n", i);
	#endif
}

/* Copy row i of B_hat matrix of size B->gSize to row vector */
void get_B_hat_row(const struct _modmat *B, num i, double *row){
	double *A_i, *K_i, *r;
	double f_i=0;
	num gSize=B->gSize;
	#ifdef DEBUG
	printf("BEGIN: get_B_hat_row=%d\n", i);
	#endif
	A_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(A_i!=NULL,MEM_ALLOC_ERROR)
	K_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(K_i!=NULL,MEM_ALLOC_ERROR)
	B->A->get_row(B->A, i, A_i);
	get_K_row(B, i, K_i);
	/* Compute entire row i of B_hat[g] and f_i */
	for (r=row; r<row+gSize; r++){
		*r=*A_i -*K_i;
		f_i += (*A_i++ - *K_i++);
	}
	*(row+i) -= f_i;
	free(A_i-gSize);
	free(K_i-gSize);
	#ifdef DEBUG
	printf("SUCCESS: get_B_hat_row=%d\n", i);
	#endif
}

double sum_of_abs(double *row, num n){
	double *p, res=0;
	for (p=row; p<row+n; p++)
		res+=fabs(*p);
	return res;
}

/* helping function, populate resK and resM with sub-vector K and M aligned with the subgroup g*/
void compute_K_and_currM(modMat *Bg, int_vector K, Subgroup g, num gSize){
    int_vector Kg = (int_vector)malloc(sizeof(num) * gSize), p, q;
    num Mg = 0;
	VERIFY(Kg!=NULL, MEM_ALLOC_ERROR)
	for (p=Bg->g, q=g; p<Bg->g+gSize; p++, q++)
		*p=*q;
	q=Bg->g;
    for (p = Kg ; p < gSize + Kg ; p++){
        *p = K[*q];
        Mg += *p;
        q++;
    }
	for (p=Bg->K, q=Kg; p<Bg->K+gSize; p++, q++)
		*p=*q;
    /*Bg->K = Kg;*/
    Bg->currM = Mg;
}

/* Constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */

modMat *create_Sub_Matrix(modMat *B, Subgroup g, num sizeG, boolean use_linked_impl){
    modMat *Bg = allocate_mod_mat(sizeG);
    VERIFY(Bg!=NULL,MEM_ALLOC_ERROR)
	if (use_linked_impl==USE_LINKED)
		Bg->A = create_sub_sparse_matrix_linked(B->A, g, sizeG, Bg->spmatSize);
	else
		Bg->A = create_sub_sparse_matrix_generic(B->A, g, sizeG, Bg->spmatSize);
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
    compute_K_and_currM(Bg, B->K, g, sizeG);
	#ifdef DEBUG
	printf("SUCCESS: create_Sub_Matrix - get K and M\n");
	#endif
	set_1_norm(Bg);
	/*Bg->one_norm=B->one_norm;*/
	#ifdef DEBUG
	printf("SUCCESS: create_Sub_Matrix\n");
	#endif
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) 
 * */
void mult_K(const modMat *Bg, const double *v, double *res){
    int_vector K = Bg->K;
    num origM = Bg->M, sizeG = Bg->gSize, *ki;
    double dot = 0;
	#ifdef DEBUG2
	printf("BEGIN: mult_K of size %d\n", sizeG);
	#endif
    VERIFY(res != NULL,NULL_POINTER_ERROR)
    for (ki = K ; ki < sizeG + K ; ki++, v++){
        dot += (*ki) * (*v);
    }
    for (ki = K ; ki < sizeG + K ; ki++, res++){
      *res = (*ki) * dot * (1 / (double)origM);
    }
	#ifdef DEBUG2
	printf("SUCCESS: mult_K of size %d\n", sizeG);
	#endif
}


/** part of the modularity matrix multiplication for power iteration, 
 * 	multiplying the 2 matrices (f_i - ||C||) * I 
 * */
void mult_F_and_C(const modMat *Bg, const double *v, double *res, boolean shift){
    int_vector K = Bg->K;
    int_vector spmatSize = Bg->spmatSize;
    num M = Bg->currM, origM = Bg->M ,sizeG = Bg->gSize, *ki;
    double fi, shiftNorm;
	#ifdef DEBUG
	printf("BEGIN: mult_F_and_C of size %d\n", sizeG);
	#endif
	shiftNorm = shift ? Bg->one_norm : 0;
    for (ki = K ; ki < sizeG + K ; ki++, v++, res++, spmatSize++){
      fi = (*spmatSize) - (((*ki) * M) / (double)origM);
      *res = (fi - shiftNorm)  * (*v);
    }
	#ifdef DEBUG
	printf("SUCCESS: mult_F_and_C of size %d\n", sizeG);
	#endif
}


/* Implements multiplication of B_hat[g] with a vector by
 * using several mult. functions and adding results together 
 * */
void mult_B_hat_g(const struct _modmat *Bg, const double *v, double *result, boolean shift){
	vector tmp1, tmp2, p;
	num gSize = Bg->gSize;
	#ifdef DEBUG
	printf("BEGIN: mult_B_hat_g of size %d\n", gSize);
	#endif
	tmp1=(vector)malloc(sizeof(double)*gSize);
	VERIFY(tmp1!=NULL,MEM_ALLOC_ERROR)
	Bg->A->mult(Bg->A,v,result);
	
	tmp1=(vector)malloc(sizeof(double)*gSize);
	VERIFY(tmp1!=NULL,MEM_ALLOC_ERROR)
	mult_K(Bg, v, tmp1);

	tmp2=(vector)malloc(sizeof(double)*gSize);
	VERIFY(tmp2!=NULL,MEM_ALLOC_ERROR)
	mult_F_and_C(Bg, v, tmp2, shift);

	for (p=result; p < result + gSize; p++)
		*p -= (*(tmp1++) + *(tmp2++));
	
	free(tmp2-gSize);
	free(tmp1-gSize);	
	#ifdef DEBUG
	printf("SUCCESS: mult_B_hat_g of size %d\n", gSize);
	#endif
}

/*
void get_B_row_generic(modMat *B, num i, double *row, boolean shift){
	vector e_i;
	#ifdef DEBUG
	printf("BEGIN: get_B_row_generic, row %d\n", i);
	#endif
	get_basis_unit_vec(&e_i, i, B->gSize);
	mult_B_hat_g(B, e_i, row, shift);
	free(e_i);
	#ifdef DEBUG
	printf("SUCCESS: get_B_row_generic, row %d\n", i);
	#endif
}
*/

void set_1_norm(modMat *B){
	/*
	 * By symmetry, row i == column i.
	 * Thus, any operation mapped (entry-wise) on a row of B is equivalent to the same operation on a column of B.
	 */
	num i, gSize=B->gSize;
	double tmp=0, max=0;
	vector B_i;
	#ifdef DEBUG
	printf("BEGIN: set_1_norm, gSize=%d\n",gSize);
	#endif
	B_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(B_i!=NULL,MEM_ALLOC_ERROR)
	for (i=0; i<gSize; i++){
		B->get_row(B,i,B_i);
		/*get_B_row_generic(B, i, B_i, NO_SHIFT);*/
		tmp=sum_of_abs(B_i, gSize);
		if (tmp>max)
			max=tmp;
	}
	free(B_i);
	#ifdef DEBUG
	printf("SUCCESS: set_1_norm\n");
	#endif
	B->one_norm=max;
}

/*allocate new ModMat of size n*/
modMat* allocate_mod_mat(num n){
	modMat *rep=(modMat*)malloc(sizeof(modMat));
	VERIFY(rep!=NULL,MEM_ALLOC_ERROR)
	rep->gSize=n;
	if (USE_LINKED)
		rep->A = spmat_allocate_list(n);
	/* Assume 2*m := count of non-zero elements required for spmat arrays impl. */ 
	/*
	else
		rep->A = spmat_allocate_array((int)n,(int)2*m);
	*/
	VERIFY(rep->A != NULL,MEM_ALLOC_ERROR)

	rep->K=(int_vector)malloc(n*sizeof(num));
	VERIFY(rep->K != NULL,MEM_ALLOC_ERROR)
	
	rep->spmatSize=(int_vector)malloc(n*sizeof(num));
	VERIFY(rep->spmatSize != NULL,MEM_ALLOC_ERROR)
	
	rep->g=(Subgroup)malloc(n*sizeof(num));
	VERIFY(rep->g != NULL,MEM_ALLOC_ERROR)
	

	rep->free=free_mod_mat;
	rep->mult=mult_B_hat_g;
	rep->get_row=get_B_hat_row;

	#ifdef DEBUG
	printf("SUCCESS: Allocated %d-sized modmat B resources\n",(int)rep->gSize);
	#endif

	return rep;
}
