#include "modMat.h"
#define DEBUG


void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->g);
	free(B);
}


void get_K_row(const modMat *B, num i, double *row){
	int_vector p;
	num k_i;
	k_i=*(B->K+i);
	for (p = B->K; p < B->K + B->gSize; p++)
		*row++=(k_i)*((double)((*p)/B->M));
}

/* Copy row i of B_hat matrix of size B->gSize to row vector */
void get_B_hat_row(const struct _modmat *B, num i, double *row){
	double *A_i, *K_i, *r;
	double f_i=0;
	num gSize=B->gSize;
	A_i=(vector)malloc(gSize*sizeof(double));
	VERIFY(A_i!=NULL,MEM_ALLOC_ERROR)
	K_i=(vector)malloc(gSize*sizeof(double));
	VERIFY (K_i!=NULL,MEM_ALLOC_ERROR)
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
	/*printf("SUCCESS: get_B_hat_row\n");*/
}

double sum_of_abs(double *row, num n){
	double *p, res=0;
	for (p=row; p<row+n; p++)
		res+=fabs(*p);
	return res;
}

/* helping function, populate resK and resM with sub-vector K and M aligned with the subgroup g*/
void gen_K_and_currM(int_vector K, Subgroup g, num gSize, modMat *Bg){
    int_vector Kg = (int_vector)malloc(sizeof(num) * gSize), p;
    num Mg = 0;
	VERIFY(Kg!=NULL, MEM_ALLOC_ERROR)
    for (p = Kg ; p < gSize + Kg ; p++){
        *p = K[*g];
        Mg += *p;
        g++;
    }
    Bg->K = Kg;
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
    Bg->g = g;
	Bg->M = B->M; /* For any submatrix created, keep the original M of the network */
    gen_K_and_currM(B->K, g, sizeG, Bg);
	#ifdef DEBUG
	printf("SUCCESS: create_Sub_Matrix get K and M\n");
	#endif
	/*set_1_norm(Bg);*/
	Bg->one_norm=B->one_norm;
	#ifdef DEBUG
	printf("SUCCESS: create_Sub_Matrix\n");
	#endif
    return Bg;
}

/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the degrees matrix (k_i * k_j / M) 
 * */
void mult_K(modMat *Bg, double *v, double *res){
    int_vector K = Bg->K;
    num origM = Bg->M, sizeG = Bg->gSize, *ki;
    double dot = 0;
    VERIFY(res != NULL,NULL_POINTER_ERROR)
    for (ki = K ; ki < sizeG + K ; ki++){
        dot += (*ki) * (*v);
        v++;
    }
	v-=sizeG;
    for (ki = K ; ki < sizeG + K ; ki++){
      *res = (*ki) * dot * (1 / (double)origM);
      res++;
    }
}


/** part of the modularity matrix multiplication for power iteration, 
 * multiplying the 2 matrices (f_i - ||C||) * I 
 * */
void mult_F_and_C(modMat *Bg, double *v, boolean shift, double *res){
    int_vector K = Bg->K;
    int_vector spmatSize = Bg->spmatSize;
    num M = Bg->currM, origM = Bg->M ,sizeG = Bg->gSize, *ki;
    double fi, shiftNorm;
	shiftNorm = shift ? Bg->one_norm : 0;
    for (ki = K ; ki < sizeG + K ; ki++){
      fi = (*spmatSize) - (((*ki) * M) / origM);
      *res = (fi - shiftNorm)  * (*v) ;
      v++;
      res++;
      spmatSize++;
    }
}


/* Implements multiplication of B_hat[g] with a vector by
 * using several mult. functions and adding results together 
 * */
void mult_B_hat_g(modMat *Bg, double *v, double *result, boolean shift){
	double *tmp1, *tmp2, *tmp3, *p;
	num gSize = Bg->gSize;
	tmp1=(double*)malloc(sizeof(double)*gSize);
	VERIFY(tmp1!=NULL,MEM_ALLOC_ERROR)
	tmp2=(double*)malloc(sizeof(double)*gSize);
	VERIFY(tmp2!=NULL,MEM_ALLOC_ERROR)
	tmp3=(double*)malloc(sizeof(double)*gSize);
	VERIFY(tmp3!=NULL,MEM_ALLOC_ERROR)

	Bg->A->mult(Bg->A,v,tmp1);
	mult_K(Bg, v, tmp2);
	mult_F_and_C(Bg, v, shift, tmp3);

	for (p=result; p < result + gSize; p++)
		*p = *(tmp1++) - *(tmp2++) - *(tmp3++);
	
	free(tmp3-gSize);
	free(tmp2-gSize);
	free(tmp1-gSize);
}


void get_B_row_generic(modMat *B, num i, double *row, boolean shift){
	vector e_i;
	#ifdef DEBUG
	printf("BEGIN: get_B_row_generic, row %d\n", i);
	#endif
	e_i=(vector)calloc(B->gSize,sizeof(double));
	VERIFY(e_i!=NULL,MEM_ALLOC_ERROR)
	*(e_i+i)=1;
	mult_B_hat_g(B, e_i, row, shift);
	free(e_i);
	#ifdef DEBUG
	printf("SUCCESS: get_B_row_generic, row %d\n", i);
	#endif
}

void set_1_norm(modMat *B){
	/*
	 * By symmetry, row i == column i.
	 * Thus, any operation mapped (entry-wise) on a row of B is equivalent to the same operation on a column of B.
	 */
	num i;
	double tmp=0, max=0;
	vector B_i;
	#ifdef DEBUG
	printf("BEGIN: set_1_norm\n");
	#endif
	B_i=(vector)malloc(B->gSize*sizeof(double));
	VERIFY(B_i!=NULL,MEM_ALLOC_ERROR)
	for (i=0; i<B->gSize; i++){
		B->get_row(B,i,B_i);
		/*get_B_row_generic(B, i, B_i, NO_SHIFT);*/
		tmp=sum_of_abs(B_i, B->gSize);
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

	rep->get_row=get_B_hat_row;

	#ifdef DEBUG
	printf("SUCCESS: Allocated %d-sized modmat B resources\n",(int)rep->gSize);
	#endif

	return rep;
}
