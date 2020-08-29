#include "spmat.h"
#define DEBUG

typedef struct linked_list {
	DATA val;
	int col;
	struct linked_list *next;
} linked_list;

typedef struct linked_list Node;

/* insert a new node after head */
void add_after(Node *head, int i){
	head->next = (Node*)malloc(sizeof(Node));
	VERIFY(head->next != NULL,MEM_ALLOC_ERROR)
	head->next->col = i;
	head->next->val = 1;
}

/* adding a row to the sub sparse matrix */
num add_row_to_sub_linked(Subgroup col, num row, Subgroup g, int i, Node *AHead, Node *subTail, int n){
	num rowSize = 0;
	while (col < n + g && AHead != NULL){
		if (row == *col || (num)AHead->col > *col){
			col++;
			i++;
		} 
		else if ((num)AHead->col == *col){
			add_after(subTail, i);
			rowSize++;
			subTail = subTail->next;
			col++;
			i++;
			AHead = AHead->next;
		}
		else{
			AHead = AHead->next;
		}	
	}
	return rowSize;
}

/* creating a sub sparse matrix for Algorithm 2 */
spmat* create_sub_sparse_matrix_linked(spmat *A, Subgroup g, int n , int_vector spmatSize){
	spmat *sub = spmat_allocate_list(n);
	Node **subSparse = sub->private, **Asparse = A->private, *AHead, *subHead, *tmp;
	Subgroup row;
	for (row = g ; row < (num)n + g ; row++){
		subHead = (Node*)malloc(sizeof(Node));
		VERIFY(subHead != NULL,MEM_ALLOC_ERROR)
		AHead = Asparse[*row];
		*spmatSize = add_row_to_sub_linked(g, *row, g, 0, AHead, subHead, n);
		tmp = subHead;
		subHead = subHead->next;
		free(tmp);
		*subSparse = subHead;
		subSparse++;
		spmatSize++;
	}
	return sub;
}

void add_row_linked(spmat *A, const double *row, int i){
	Node *head = NULL, *tail;
	int j = 0, n = A->n;
	const DATA *rowPtr = row;
	while (j < n && *rowPtr == 0){ /* getting to the first non-zero value and adding it to the first node*/
		rowPtr++;
		j++;
	}
	if (j < n){
		head = (Node*) malloc (sizeof(Node));
		VERIFY(head != NULL,MEM_ALLOC_ERROR)
		head->val = *rowPtr;
		head->col = j;
		tail = head;
		rowPtr++;
		j++;
		while (j < n){ /* add to tail*/
			while (j < n && *rowPtr == 0){
				rowPtr++;
				j++;
			}
			if (j < n){
				tail->next = (Node*) malloc (sizeof(Node));
				VERIFY(tail->next != NULL,MEM_ALLOC_ERROR)
				tail = tail->next;
				tail->val = *rowPtr;
				tail->col = j;
				rowPtr++;
				j++;
			}
		}
		tail->next = NULL; /* end of list*/
	}
	((Node**)(A->private))[i] = head;
}
/* helping function, dot product of a row of the sparse matrix and the input vector*/
DATA dot_product(Node* row, const DATA *col){
	Node* head = row;
	DATA sum = 0;
	while(head != NULL){
		sum += (head->val) * (col[head->col]);
		head = head->next;
	}
	return sum;
}

void mult_linked(const spmat *A, const DATA *v, DATA *result){
	Node **mat = (Node**)(A->private);
	DATA *resPtr = result;
	int n = A->n, i;
	for (i = 0; i < n ; i++){
		if(*mat == NULL){
			*resPtr = 0; /* if there is a row with zero values only*/
		}
		else{
			*resPtr = dot_product(*mat,v); /* dot product of each row of the sparse matrix and the input vector*/
		}
		resPtr++;
		mat++;
	}
}


void delete_list(Node *head){
	if (head != NULL){
		delete_list(head->next);
		free(head);
	}
}

void free_linked(spmat *A){
	Node** mat = (Node**)(A->private);
	int n = A->n, i;
	for (i = 0 ; i < n ; i++){
		delete_list(*mat);
		mat++;
	}
	free(A->private);
	free(A);
}



spmat* spmat_allocate_list(int n){
	spmat *A;
	A = (spmat*) malloc (sizeof(spmat));
	VERIFY(A != NULL,MEM_ALLOC_ERROR)
	A->n = n;
	A->add_row = add_row_linked;
	A->free = free_linked;
	A->mult = mult_linked;
	/*A->get_row = get_spmat_row_linked;*/
	A->private =  (Node**) malloc (sizeof(Node*) * (A->n));
	VERIFY(A->private != NULL,MEM_ALLOC_ERROR)
	return A;
}

/*** ARRAY MATRIX IMPLEMENTATION - TO BE DELETED OR MODIFIED ***/

typedef struct _arrmat
{
	scalar* values;
	int* colind;
	int* rowptr;
} arraymat;

void add_row_arrays(struct _spmat *A, const double *row, int i)
{
	arraymat *imp=((arraymat*)(A->private));
	scalar *valOffset=imp->values;
	const double *vecEntry;
	int *colOffset=imp->colind, *rp=imp->rowptr+i, rowNNZ=0;
	VERIFY(i < A->n,OUT_OF_BOUNDS_ERROR)

	/* Find the last position in values array that stores a valid non-zero value.
	*/
	valOffset += i>0 ? *rp : 0;
	colOffset += valOffset-(imp->values);

	/* Add non-zero entries of row vector into A's values and column ind. arrays.
	 */
	for (vecEntry=row; vecEntry<row+A->n; vecEntry++)
		if (*vecEntry!=0.0){
			rowNNZ++;
			*valOffset++=*vecEntry;
			*colOffset++=vecEntry-row;
		}

	/* A nonzero difference between consecutive row pointers indicates that current row r is populated with
	 * non-zero elements. Otherwise, the first instance of a non-zero value applies only to a later row.
	 * The difference is also the actual number of non-zero values in that row.
	 */
	*(rp+1)=*rp+rowNNZ;
}


void get_row_arrays(const struct _spmat *A, int i, double *result){
	arraymat *imp=((arraymat*)(A->private));
	double *p, *v=imp->values+*(imp->rowptr+i);
	int *src_col = imp->colind+*(imp->rowptr+i), dest_col;
	for (p=result; p<result+A->n; p++){
		dest_col=p-result;
		if (dest_col==*src_col){
			*p=*v;
			v++;
			src_col++;
		}
		else
			*p=0;
	}
}

void free_arrays(struct _spmat *A)
{
	arraymat *imp=((arraymat*)(A->private));
	free(imp->values);
	free(imp->colind);
	free(imp->rowptr);
	free(A->private);
	free(A);
}

void mult_arrays(const struct _spmat *A, const double *v, double *result)
{
	arraymat *imp=((arraymat*)(A->private));
	int *r, *j, rowNNZ;
	scalar sum, *u;
	for (r=imp->rowptr; r<imp->rowptr+A->n; r++)
	{
		sum=0.0;
		rowNNZ=*(r+1)-*r;
		u=imp->values+(*r);
		j=imp->colind+(*r);
		while (rowNNZ-->0)
			sum+=(*u++)*(*(v+*j++));
		*result++=sum;
	}
}

/** TO BE TESTED.
 * 
 * 	Generic get_row for spmat, independently of implementation (lists or arrays).
 * 	Utilizes A->mult with a standard vector (v[j]==1 if i=j, else v[j]==0).	
 *  
 */
/*
void get_spmat_row_with_mult(const struct _spmat *A, double *row, int i){
	double *basis;
	basis=(double*)calloc(A->n,sizeof(double));
	VERIFY(basis!=NULL,MEM_ALLOC_ERROR)
	*(basis+i)=1;
	A->mult(A,basis,row);
	free(basis);
}
*/

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz)
{
	spmat *ret=(spmat*)malloc(sizeof(spmat));
	arraymat *imp=(arraymat*)malloc(sizeof(arraymat));
	VERIFY(ret != NULL,MEM_ALLOC_ERROR)
	VERIFY(imp == NULL,MEM_ALLOC_ERROR)
	imp->values = (scalar*)malloc(nnz*sizeof(scalar));
	VERIFY(imp->values != NULL,MEM_ALLOC_ERROR)
	imp->colind = (int*)malloc(nnz*sizeof(int));
	VERIFY(imp->colind != NULL,MEM_ALLOC_ERROR)
	imp->rowptr = (int*)malloc((n+1)*sizeof(int));
	VERIFY(imp->rowptr != NULL,MEM_ALLOC_ERROR)
	ret->n=n;
	ret->get_row=get_row_arrays;
	ret->add_row=add_row_arrays;
	ret->free=free_arrays;
	ret->mult=mult_arrays;
	ret->private=imp; /* pointer assignment */
	return ret;
}

num get_sub_matrix_nnz_arrays(spmat *A, Subgroup g, int sizeG){
	arraymat *imp=((arraymat*)(A->private));
	int *cols=imp->colind, *t;
	Subgroup p;
	num cnt=0;
	for (p=g; p<g+sizeG; p++){
		for (t=cols; t<cols+A->n; t++)
			cnt += ((int)*p)==*t ? 1 : 0;
	}
	return cnt;
}

/**
 * Can be used as a generic interface function of spmat module (for both impl.)
 */
void get_sub_row_arrays(spmat *A, int i, double *dest, Subgroup g, num sizeG, int_vector spmatSize){
	double *v, sub_ij;
	Subgroup q;
	v=(double*)malloc(A->n*sizeof(double));
	VERIFY(v!=NULL, MEM_ALLOC_ERROR)
	A->get_row(A,i,v);
	for (q=g; q<g+sizeG; q++){
		sub_ij=*(v+*q);
		*(dest++)=sub_ij;
		*(spmatSize+i)+=(num)sub_ij;
	}
	free(v);
}
/* Creating a sub sparse matrix for Algorithm 2.
 * Can be used as a generic interface function of spmat module (for both impl.)
 * 
 * If impl_flag==1, uses linked-list implementation. Otherwise, use arrays impl.
 */
spmat* create_sub_sparse_matrix_array(spmat *A, Subgroup g, int sizeG , int_vector spmatSize /*, int impl_flag*/){
	spmat *sub;
	double *sub_row;
	Subgroup p;
	/*if (impl_flag)
		sub = spmat_allocate_linked(sizeG);
	else */
		sub = spmat_allocate_array(sizeG,get_sub_matrix_nnz_arrays(A, g, sizeG));
	sub_row=(double*)malloc(sizeG*sizeof(double));
	VERIFY(sub_row!=NULL, MEM_ALLOC_ERROR)
	for (p=g; p<g+sizeG; p++){
		get_sub_row_arrays(A,*p,sub_row,g,sizeG,spmatSize);
		sub->add_row(sub,sub_row,p-g);
	}
	free(sub_row);
	return sub;
}

