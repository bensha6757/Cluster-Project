#include "modMat.h"


void free_mod_mat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->g);
	free(B);
}


void get_adj_row(const modMat *B, size_t i, double *row){
	B->A->get_row(B->A, i, row);
}

void get_K_row(const modMat *B, size_t i, double *row){
	int_vector p;
	size_t k_i;
	k_i=*(B->K+i);
	for (p = B->K; p < B->K + B->gSize; p++)
		*row++=(k_i)*((double)((*p)/B->M));
}

/* Copy row i of B_hat matrix of size B->gSize to row vector */
void get_B_hat_row(const struct _modmat *B, size_t i, double *row){
	double *A_i, *K_i;
	double f_i=0;
	size_t j;
	VERIFY(i<B->gSize,OUT_OF_BOUNDS_ERROR)
	A_i=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(A_i!=NULL,MEM_ALLOC_ERROR)
	K_i=(double*)malloc(B->gSize*sizeof(double));
	VERIFY (A_i!=NULL,MEM_ALLOC_ERROR)
	get_adj_row(B,i,A_i);
	get_K_row(B,i,K_i);
	/* Compute entire row of B_hat[g] and f_i */
	for (j=0; j<B->gSize; j++){
		*row++=*A_i++ -*K_i++;
		f_i += (*A_i++ - *K_i++);
	}
	row-=B->gSize;
	*(row+i)-=f_i;
	A_i-=B->gSize;
	free(A_i);
	K_i-=B->gSize;
	free(K_i);
}

double sum_of_abs(double *row, size_t n){
	double *p, res=0;
	for (p=row; p<row+n; p++)
		res+=fabs(*p);
	return res;
}


void set_1_norm(modMat *B){
	/*
	 * By symmetry, row i == column i.
	 * Thus, any operation mapped on a row of B is equivalent to the same operation on columns.
	 */
	size_t i;
	double tmp=0, max=0;
	double *B_i=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(B_i!=NULL,MEM_ALLOC_ERROR)
	for (i=0; i<B->gSize; i++){
		B->get_row(B,i,B_i);
		tmp=sum_of_abs(B_i, B->gSize);
		if (tmp>max)
			max=tmp;
	}
	free(B_i);
	B->one_norm=max;
}

/* helping function, populate resK and resM with sub-vector K and M aligned with the subgroup g*/
void genKandM(int_vector K, Subgroup g, size_t gSize, modMat *Bg){
    int_vector Kg = (int_vector)malloc(sizeof(size_t) * gSize), p;
    size_t Mg = 0;
	VERIFY(Kg!=NULL,MEM_ALLOC_ERROR)
    for (p = Kg ; p < gSize + Kg ; p++){
        *p = K[*g];
        Mg += *p;
        g++;
    }
    Bg->K = Kg;
    Bg->M = Mg;
}

/* helping function, generates M, which is the sum of all degrees in K*/
size_t genM(int_vector K, size_t sizeG){
    size_t M = 0, *p;
    for (p = K ; p < sizeG + K ; p++){
      M += *p;
    }
    return M;
}

/* constructor, creating a new sub matrix B_hat[g]. 
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */

modMat *create_Sub_Matrix(modMat *B, Subgroup g, size_t sizeG, boolean impl_flag){
    modMat *Bg = allocate_mod_mat(sizeG);
    VERIFY(Bg!=NULL,MEM_ALLOC_ERROR)
	if (impl_flag==USE_LINKED)
    	Bg->A = create_sub_sparse_matrix_linked(B->A, g, sizeG, Bg->spmatSize);
	else
		Bg->A = create_sub_sparse_matrix_array(B->A, g, sizeG, Bg->spmatSize);
    Bg->g = g;
    genKandM(B->K, g, sizeG, Bg);
	set_1_norm(Bg);
    return Bg;
}

/* part of the modularity matrix multiplication for power iteration, multiplying the degrees matrix (k_i * k_j / M) */
void mult_K(modMat *B, modMat *Bg, double *v, double *res){
    int_vector K = Bg->K;
    size_t origM = B->M, sizeG = Bg->gSize, *ki;
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


/* part of the modularity matrix multiplication for power iteration, multiplying the 2 matrices (f_i - ||C||) * I */
void mult_F_and_C(modMat *B, modMat *Bg, double *v, boolean shift, double *res){
    int_vector K = Bg->K;
    int_vector spmatSize = Bg->spmatSize;
    size_t M = Bg->M, origM = B->M ,sizeG = Bg->gSize, *ki;
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
 * using several mult. functions and adding results together */
void mult_B_hat_g(modMat *B, modMat *Bg, double *v, double *result){
	double *tmp1, *tmp2, *tmp3, *p;
	tmp1=(double*)malloc(sizeof(double)*B->gSize);
	VERIFY(tmp1!=NULL,MEM_ALLOC_ERROR)
	tmp2=(double*)malloc(sizeof(double)*B->gSize);
	VERIFY(tmp2!=NULL,MEM_ALLOC_ERROR)
	tmp3=(double*)malloc(sizeof(double)*B->gSize);
	VERIFY(tmp3!=NULL,MEM_ALLOC_ERROR)
	Bg->A->mult(Bg->A,v,tmp1);
	mult_K(B, Bg, v, tmp2);
	mult_F_and_C(B, Bg, v, (Bg->gSize != B->gSize), tmp3);
	for (p=result; p<result+Bg->gSize; p++)
		*p=*tmp1++ - *tmp2++ - *tmp3++;
	tmp3-=Bg->gSize;
	free(tmp3);
	tmp2-=Bg->gSize;
	free(tmp2);
	tmp1-=Bg->gSize;
	free(tmp1);
}

/*allocate new ModMat of size n*/
modMat* allocate_mod_mat(int n){
	modMat *rep=(modMat*)malloc(sizeof(modMat));
	VERIFY(rep!=NULL,MEM_ALLOC_ERROR)
	rep->gSize=n;

	if (USE_LINKED)
		rep->A = spmat_allocate_list(n);
	VERIFY(rep->A != NULL,MEM_ALLOC_ERROR)

	rep->K=(int_vector)malloc(n*sizeof(size_t));
	VERIFY(rep->K != NULL,MEM_ALLOC_ERROR)

	rep->g=(Subgroup)malloc(n*sizeof(size_t));
	VERIFY(rep->g != NULL,MEM_ALLOC_ERROR)
	
	rep->free=free_mod_mat;

	rep->get_row=get_B_hat_row;

	return rep;
}


/* 	DEPRECATED
 *
 *	Multiply submatrix B[g] with vector u;
 *  This is an implementation for mult method in ModMat type struct.
 *	For original B of the whole network, sumK == B->M 
 *  
 *  spmat's mult must be modified to support multiplying B[g]
 */

/*
void multSubgroupMatrix(modMat *Bg, Subgroup g, double *u, double *res){
	double *resAu, *resKu, sumK=0;
	vector j, *i, t, *tempg=g;
	resAu=(double*)malloc(Bg->n*sizeof(double));
	Bg->A->mult(Bg->A,u,resAu); 
	resKu=(double*)malloc(Bg->n*sizeof(double));
	for (j=Bg->K; j<Bg->K+Bg->n; j++)
		if (*tempg==j-Bg->K && *tempg<Bg->n){
			sumK += (*j);
			tempg++;
		}
	
	tempg=g;
	for (i=Bg->K; i<Bg->K+Bg->n; i++)
		if (*tempg==i-Bg->K && *tempg<Bg->n){
			*resKu++ = ((double)((*i)*sumK) / Bg->M)*(*u++);
			tempg++;
		}
	resKu-=Bg->n;
	for (t=0; t<Bg->n; t++)
		*res++=(*resAu++)-(*resKu++);
	free(resAu);
	free(resKu);
}
*/