#include "spmat.h"

typedef struct linked_list {
	DATA val;
	int col;
	struct linked_list *next;
} linked_list;

typedef struct linked_list Node;


 /************ DEPRECATED ************/

 /*
void get_basis_unit_vec(vector *e_i, num i, num n){
	*e_i=(vector)calloc(n,sizeof(double));
	VERIFY(e_i!=NULL,MEM_ALLOC_ERROR)
	*(*e_i+i)=1;
}
*/


/** GENERIC FUNCTION, INDEP. OF SPMAT IMPL.
 *  Fetch row i of spmat A, by utilizing A->mult with a standard vector (v[j]==1 if i=j, else v[j]==0)
 */

/*
void get_spmat_row_generic(const struct _spmat *A, int i, double *row){
	vector e_i;
	#ifdef DEBUG_SPMAT
	printf("BEGIN: get_spmat_row_generic=%d\n", i);
	#endif
	get_basis_unit_vec(&e_i, i, A->n);
	A->mult(A,e_i,row);
	free(e_i);
	#ifdef DEBUG_SPMAT
	printf("SUCCESS: get_spmat_row_generic=%d\n", i);
	#endif
}
*/

/*** LINKED LIST MATRIX IMPLEMENTATION ***/

void get_row_linked(const struct _spmat *A, int i, vector row){
	Node **mat = (Node**)(A->private), *head;
	int n = A->n, j;
	head = mat[i];
	for (j = 0; j < n && head != NULL; row++ ,j++){
		if (j < head->col){
			*row = 0;
		}
		else{
			*row = 1;
			head = head->next;
		}
	}
	while (j < n)
	{
		*(row++) = 0;
		j++;
	}
}


/* insert a new node after head */
void add_after(Node *head, int i){
	head->next = (Node*)malloc(sizeof(Node));
	VERIFY(head->next != NULL,MEM_ALLOC_ERROR)
	head->next->col = i;
	head->next->val = 1;
}


/* adding a row to the sub sparse matrix */
num add_row_to_sub_linked(Subgroup col, num row, Subgroup g, Node *AHead, Node *subHead, int n){
	num rowSize = 0, i=0;
	while (col < n + g && AHead != NULL){
		if (row == *col || (num)AHead->col > *col){
			col++;
			i++;
		} 
		else if ((num)AHead->col == *col){
			add_after(subHead, i);
			rowSize++;
			subHead = subHead->next;
			col++;
			i++;
			AHead = AHead->next;
		}
		else{
			AHead = AHead->next;
		}	
	}
	subHead->next = NULL; /* end of list */
	return rowSize;
}


/* creating a sub sparse matrix for Algorithm 2 */
spmat* create_sub_sparse_matrix_linked(spmat *A, Subgroup g, int n , int_vector spmatSize){
	spmat *sub = spmat_allocate_list(n);
	Node **subMat = sub->private, **AMat = A->private, *AHead, *subHead, *tmp;
	Subgroup row;
	#ifdef DEBUG_SPMAT
	printf("BEGIN: create_sub_sparse_matrix_linked of size %d\n",n);
	#endif
	for (row = g ; row < g + n ; row++){
		subHead = (Node*)malloc(sizeof(Node));
		VERIFY(subHead != NULL,MEM_ALLOC_ERROR)
		tmp = subHead;
		AHead = AMat[*row];
		if (AHead != NULL){
			*spmatSize = add_row_to_sub_linked(g, *row, g, AHead, subHead, n);
		}
		else{
			*spmatSize = 0;
			subHead->next = NULL;
		}
		spmatSize++;
		subHead = subHead->next;
		free(tmp);
		*(subMat++) = subHead;
	}
	#ifdef DEBUG_SPMAT
	printf("SUCCESS: create_sub_sparse_matrix_linked of size %d\n",n);
	#endif
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
DATA dot_product_linked(Node* row, const DATA *col){
	Node* head = row;
	DATA sum = 0;
	while(head != NULL){
		sum += (head->val) * (col[head->col]);
		head = head->next;
	}
	return sum;
}

void mult_linked(const struct _spmat *A, const double *v, double *result){	
	Node **mat = (Node**)(A->private);
	DATA *resPtr = result;
	int n = A->n, i;
	#ifdef DEBUG_SPMAT
	printf("BEGIN: mult_linked of size %d\n",n);
	#endif
	for (i = 0; i < n ; i++){
		if (*mat == NULL){
			*resPtr = 0; /* if there is a row with zero values only*/
		}
		else {
			*resPtr = dot_product_linked(*mat,v); /* dot product of each row of the sparse matrix and the input vector*/
		}
		resPtr++;
		mat++;
	}
	#ifdef DEBUG_SPMAT
	printf("SUCCESS: mult_linked of size %d\n",n);
	#endif
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
	A->get_row = get_row_linked;
	A->private =  (Node**) malloc (sizeof(Node*) * (A->n));
	VERIFY(A->private != NULL,MEM_ALLOC_ERROR)
	return A;
}

/*** ARRAY MATRIX IMPLEMENTATION ***/

typedef struct _arrmat
{
	scalar* values;
	int_vector colind;
	int_vector rowptr;
} arraymat;

void add_row_arrays(struct _spmat *A, const double *row, int i)
{
	arraymat *imp=((arraymat*)(A->private));
	scalar *valOffset=imp->values;
	const double *vecEntry;
	num *colOffset=imp->colind, *rp=imp->rowptr+i, rowNNZ=0;
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

void get_row_arrays(const struct _spmat *A, int i, double *row){
	arraymat *imp=((arraymat*)(A->private));
	double *p, *v=imp->values+*(imp->rowptr+i);
	num *src_col = imp->colind+*(imp->rowptr+i), dest_col;
	for (p=row; p<row+A->n; p++){
		dest_col=p-row;
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
	num *r, *j, rowNNZ;
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

/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz)
{
	spmat *ret=(spmat*)malloc(sizeof(spmat));
	arraymat *imp=(arraymat*)malloc(sizeof(arraymat));
	VERIFY(ret != NULL,MEM_ALLOC_ERROR)
	VERIFY(imp != NULL,MEM_ALLOC_ERROR)
	imp->values = (scalar*)malloc(nnz*sizeof(scalar));
	VERIFY(imp->values != NULL,MEM_ALLOC_ERROR)
	imp->colind = (int_vector)malloc(nnz*sizeof(num));
	VERIFY(imp->colind != NULL,MEM_ALLOC_ERROR)
	imp->rowptr = (int_vector)malloc((n+1)*sizeof(num));
	VERIFY(imp->rowptr != NULL,MEM_ALLOC_ERROR)
	ret->n=n;
	ret->get_row=get_row_arrays;
	ret->add_row=add_row_arrays;
	ret->free=free_arrays;
	ret->mult=mult_arrays;
	ret->private=imp; /* pointer assignment */
	return ret;
}

num get_g_row_nnz_arrays(arraymat *orig, int i, Subgroup g, int sizeG){
	num *r = orig->rowptr+i, *cols=orig->colind, *orig_c=cols+*r;
	num cnt=0, row_nnz=*(r+1)-*r;
	Subgroup p=g;
	while (orig_c < cols + *r + row_nnz && p<g+sizeG)
		if 		(*p < *orig_c)
			p++;
		else if (*p > *orig_c)
			orig_c++;
		else {
			cnt++;
			p++;
			orig_c++;
	}
	return cnt;
}

num get_g_nnz_arrays(spmat *A, Subgroup g, int sizeG, int_vector spmatSize){
	arraymat *orig=((arraymat*)(A->private));
	Subgroup p;
	num cnt=0, k;
	for (p=g; p<g+sizeG; p++){
		k=get_g_row_nnz_arrays(orig,*p,g,sizeG);
		spmatSize[p-g]=k;
		cnt+=k;
	}
		
	return cnt;
} 

void add_row_to_sub_arrays(arraymat *orig, int orig_i, arraymat *dest, int sub_i, Subgroup g, num sizeG){
	num 	k=orig->rowptr[orig_i], orig_nnz=orig->rowptr[orig_i+1]-k;
	num 	*c_orig=orig->colind+k;
	double 	*v_orig=orig->values+k;
	num 	*r_dest=dest->rowptr+sub_i;
	num		*c_dest=dest->colind, *c_p;
	double 	*v_dest=dest->values;
	Subgroup p=g;
	#ifdef DEBUG_SPMAT
	printf("BEGIN: add_row_to_sub_arrays orig %d to %d\n",orig_i,sub_i);
	#endif

	if (sub_i==0)
		*r_dest=0;
	else {
		c_dest += *r_dest;
		v_dest += *r_dest;
	}

	c_p=c_dest;
	while(c_orig < orig->colind+k+orig_nnz && p < g+sizeG){
		if (*c_orig < *p){
			c_orig++;
			v_orig++;
		}
		else if (*c_orig > *p){
			p++;
		}
		else {
			*(v_dest++)=*v_orig++;
			*(c_p++)=p-g;
			c_orig++;
			p++;
		}
	}

	*(r_dest+1)= *r_dest + (c_p-c_dest); 
	#ifdef DEBUG_SPMAT
	printf("SUCCESS: add_row_to_sub_arrays orig %d to %d\n",orig_i,sub_i);
	#endif
}


spmat* create_sub_sparse_matrix_array(spmat *A, Subgroup g, int sizeG, int_vector spmatSize){
	arraymat *orig=(arraymat*)(A->private), *imp;
	spmat *sub;
	Subgroup p;
	num sub_nnz = get_g_nnz_arrays(A, g, sizeG, spmatSize);
	#ifdef DEBUG_SPMAT
	printf("sub_nnz = %d\n", sub_nnz);
	#endif
	sub=spmat_allocate_array(sizeG, sub_nnz);
	imp=(arraymat*)(sub->private);
	for (p=g; p<g+sizeG; p++)
		add_row_to_sub_arrays(orig, *p, imp, p-g, g, sizeG);
	#ifdef DEBUG_SPMAT
	printf("SUCCESS: create_sub_sparse_matrix_array of size %d\n",sizeG);
	#endif
	return sub;
}

/************ DEPRECATED *************/

/**
 * GENERIC FUNCTION, INDEP. OF SPMAT IMPL.
 * Fetch a reduced row i of A according to subgroup g and store it in dest.
 */
void get_sub_row_generic(spmat *A, int i, double *dest, Subgroup g, num sizeG, int_vector spmatSize){
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

/** GENERIC FUNCTION, INDEP. OF SPMAT IMPL.
 *  Create a sub sparse matrix for Algorithm 2.
 */	
spmat* create_sub_sparse_matrix_generic(spmat *A, Subgroup g, int n , int_vector spmatSize){
	spmat *sub;
	double *sub_row;
	Subgroup p;
	sub = spmat_allocate_list(n);
	sub_row=(double*)malloc(n*sizeof(double));
	VERIFY(sub_row!=NULL, MEM_ALLOC_ERROR)
	for (p=g; p<g+n; p++){
		get_sub_row_generic(A,*p,sub_row,g,n,spmatSize);
		sub->add_row(sub,sub_row,p-g);
	}
	free(sub_row);
	return sub;
}
