#include "spmat.h"
#include "modMat.h"
#include "io_mem_errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

typedef size_t* vector;
typedef vector Subgroup;


void freeModMat(modMat *B){
	B->A->free(B->A);
	free(B->K);
	free(B->g);
	free(B);
}

/*allocate new ModMat of size n*/
modMat *allocateModMat(int n){
	char z='z';
	modMat *rep=(modMat*)malloc(sizeof(modMat));
	if (rep==NULL)
		exit(MEM_ALLOC_ERROR);

	rep->gSize=n;

	rep->A = spmat_allocate_list(n);
	if (rep->A==NULL)
		z='a';

	rep->K=(double*)malloc(n*sizeof(double));
	if (rep->K==NULL)
		z='b';

	rep->g=(Subgroup)malloc(n*sizeof(size_t));
	if (rep->g==NULL)
		z='c';

	switch(z){
	case 'c':	free(rep->K);
	case 'b':	rep->A->free(rep->A);
	case 'a':	free(rep);
	default : 	exit(MEM_ALLOC_ERROR);
	}
	rep->free=freeModMat;

	rep->mult=multB_hat_g;

	rep->get_row=getBhatRow;

	return rep;
}


void getAdjRow(modMat *B, size_t i, double *row){
	B->A->get_row(B->A, row, i);
}

void getKRow(modMat *B, size_t i, double *row){
	vector p;
	size_t k_i;
	k_i=*(B->K+i);
	for (p=B->K; p<B->K+B->gSize; p++)
		*row++=(k_i)*((double)((*p)/B->M));
}

void getBhatRow(modMat *B, size_t i, double *row){
	double *A_i, *K_i;
	double f_i=0;
	size_t j;
	verify(i<B->gSize);
	A_i=(double*)malloc(B->gSize*sizeof(double));
	if (A_i==NULL)
		exit(MEM_ALLOC_ERROR);
	K_i=(vector)malloc(B->gSize*sizeof(size_t));
	if (A_i==NULL){
		free(A_i);
		exit(MEM_ALLOC_ERROR);
	}
	getAdjRow(B,i,A_i);
	getKRow(B,i,K_i);
	/*Compute f_i*/
	for (j=0; j<B->gSize; j++)
		f_i += *A_i++ + *K_i++;
	for (j=0; j<B->gSize; j++)
		*row++=*A_i++ -*K_i++;
	*(row+i)-=f_i;
	free(A_i);
	free(K_i);
}

double sumOfAbs(double *row, size_t n){
	double *p, res=0;
	for (p=row; p<row+n; p++)
		res+=fabs(*p);
	return res;
}



double getOneNorm(modMat *B){
	/*
	 * By symmetry, row i == column i.
	 * Thus, any operation mapped on a row of B is equivalent to the same operation on columns.
	 */
	size_t i;
	double tmp=0, max=0;
	double B_i=(double*)malloc(B->gSize*sizeof(double));
	if (B_i==NULL)
		exit(MEM_ALLOC_ERROR);
	for (i=0; i<B->gSize; i++){
		getBhatRow(B,i,B_i);
		tmp=sumOfAbs(B_i, B->gSize);
		if (tmp>max)
			max=tmp;
	}
	free(B_i);
	return max;
}

/*********************************
 * INPUT FILE TO MEMORY FUNCTIONS*
 *********************************/

/*read bytes from file or exit if failed. used to replace assert functionality*/
void readFromFile(int *dest,unsigned int num,unsigned int size,FILE* file){
	if (fread(dest,num,size,file)!=size)
		exit(FILE_READ_ERROR);
}

/*	Convert node indices (like input file format) to a {0,1} vector of length n
 * 	that can be added to a spmat via add_row. */
void convertAdjList(size_t k, size_t n, vector adj, vector res){
	int *i, temp;
	for (i=adj; i<adj+k; i++){
		temp=*i;
		while (temp-->0)
			*res++=0;
		*res++=1;
	}
	for(temp=n-*i; temp>0; temp--)
		*res++=0;
}


void readVnumFromFile(FILE *input, size_t *n){
	readFromFile(n,sizeof(size_t),1,input);
	rewind(input);
}

void loadModMatrixFromFile(FILE *input, modMat *B){
	size_t i, currDeg;
	vector matLine, *inputNeighbours, *k, *g;
	readFromFile(i,sizeof(size_t),1,input); /*Assuming file rewinded, skip |V| */
	k=B->K;
	g=B->g;
	/* Populate B with A, K matrices, and compute M */
	for (i=0; i<B->gSize; i++){
		readFromFile(currDeg,sizeof(size_t),1,input);
		*k++=currDeg;
		B->M+=currDeg;
		inputNeighbours=(vector)malloc(currDeg*sizeof(size_t));
		matLine=(vector)malloc(B->gSize*sizeof(size_t));
		if (inputNeighbours==NULL)
			exit(MEM_ALLOC_ERROR);
		readFromFile(inputNeighbours,sizeof(size_t),currDeg,input);
		convertAdjList(currDeg, B->gSize, inputNeighbours, matLine);
		B->A->add_row(B->A,matLine,i);
		*g++=1;
		free(matLine);
		free(inputNeighbours);
	}
	rewind(input);
}

/* helping function, populate resK and resM with sub-vector K and M aligned with the subgroup g*/
void genKandM(vector K, Subgroup g, size_t gSize, modMat *Bg){
    vector Kg = (vector)malloc(sizeof(size_t) * gSize);
    size_t Mg = 0;
    for (vector p = Kg ; p-Kg < gSize ; p++){
        *p = K[*g];
        Mg += *p;
        g++;
    }
    Bg->K = Kg;
    Bg->M = Mg;
}

/* helping function, generates M, which is the sum of all degrees in K*/
size_t genM(vector K, size_t sizeG){
    size_t M = 0;
    for (size_t *p = K ; p-K < sizeG ; p++){
      M += *p;
    }
    return M;
}

/* constructor, creating a new sub matrix B_hat[g] */
modMat* create_Sub_Matrix(modMat *B, Subgroup g, size_t sizeG){
    modMat* Bg = allocateModMat(sizeG);
    verify(Bg != NULL);
    Bg->A = create_sub_sparse_matrix_linked(B->A, g, sizeG, Bg->spmatSize);
    Bg->g = g;
    genKandM(B->K, g, sizeG, Bg);
    return Bg;
}

double dot_product(vector K, double *v, size_t sizeG){
    double res = 0;
    for (size_t *ki = K ; ki-K < sizeG ; ki++){
        res += (*ki) * (*v);
        v++;
    }
    return res;
}

/* part of the modularity matrix multiplication for power iteration, multiplying the degrees matrix (k_i * k_j / M) */
void mult_K(modMat *B, modMat *Bg, double *v, double *res){
    vector K = Bg->K;
    size_t origM = B->M, sizeG = Bg->gSize;
    double dot = dot_product(K,v,sizeG);
    verify(res != NULL);
    for (size_t *ki = K ; ki-K < sizeG ; ki++){
      *res = (*ki) * dot * (1 / (double)origM);
      res++;
    }
}


/* part of the modularity matrix multiplication for power iteration, multiplying the 2 matrices (f_i - ||C||) * I */
void mult_F_and_C(modMat *B, modMat *Bg, double *v, bool shift, double *res){
    vector K = Bg->K;
    vector spmatSize = Bg->spmatSize;
    size_t M = Bg->M, origM = B->M ,sizeG = Bg->gSize;
    double fi, shiftNorm;
    if (!shift)
      shiftNorm = getOneNorm(Bg);
    for (size_t *ki = K ; ki-K < sizeG ; ki++){
      fi = (*spmatSize) - (((*ki) * M) / origM);
      *res = (fi - shiftNorm)  * (*v) ;
      v++;
      res++;
      spmatSize++;
    }
}


/*Implements multiplication of B_hat[g] with a vector by
 * using several mult. functions and adding results together */
void multB_hat_g(modMat *B, modMat Bg, const double *v, double *result){
	double *tmp1, *tmp2, *tmp3, shiftNorm, *p;
	tmp1=(double*)malloc(sizeof(double)*B->gSize);
	if (tmp1==NULL)
		exit(MEM_ALLOC_ERROR);
	tmp2=(double*)malloc(sizeof(double)*B->gSize);
	if (tmp2==NULL){
		free(tmp1);
		exit(MEM_ALLOC_ERROR);
	}
	tmp3=(double*)malloc(sizeof(double)*B->gSize);
	if (tmp3==NULL){
		free(tmp2);
		free(tmp1);
		exit(MEM_ALLOC_ERROR);
	}
	verify(tmp1!=NULL && tmp2!=NULL && tmp3!=NULL);
	Bg->A->mult(Bg->A,v,tmp1);
	mult_K(B, Bg, v, tmp2);
	mult_F_and_C(B, Bg, v, Bg->gSize == B->gSize, tmp3);
	for (p=result; p<result+Bg->gSize; p++)
		*p=*tmp1++ - *tmp2++ - *tmp3++;
	free(tmp3);
	free(tmp2);
	free(tmp1);
}


/*
 *
 *Multiply submatrix B[g] with vector u;
  This is an implementation for mult method in ModMat type struct.

void multSubgroupMatrix(modMat *Bg, Subgroup g, double *u, double *res){
	double *resAu, *resKu, sumK=0;
	vector j, *i, t, *tempg=g;

	resAu=(double*)malloc(Bg->n*sizeof(double));
	Bg->A->mult(Bg->A,u,resAu); /* spmat's mult must be modified to support multiplying B[g]

	resKu=(double*)malloc(Bg->n*sizeof(double));
	for (j=Bg->K; j<Bg->K+Bg->n; j++)
		if (*tempg==j-Bg->K && *tempg<Bg->n){
			sumK += (*j);
			tempg++;
		}
	/* For original B of the whole network, sumK == B->M

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
