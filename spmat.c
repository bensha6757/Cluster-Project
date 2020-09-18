#include "Spmat.h"

typedef struct linked_list {
	DATA val;
	int col;
	struct linked_list *next;
} linked_list;

typedef struct linked_list Node;


/*** LINKED LIST MATRIX IMPLEMENTATION ***/

double get_modularity_score_linked(const struct _spmat *A, vector s, int moved_v, int_vector K, num M){
	Node **mat = (Node**)(A->private), *head;
	num gSize = A->n;
	double dQ;
	vector d_j;
	int_vector K_j;
	num k_i = K[moved_v];
	int j, d_i;
	double k_i_M = k_i / (double)M, sum = 0;
	
	s[moved_v] *= -1;
	d_i = s[moved_v];

	head = mat[moved_v];
	for (d_j = s, K_j = K, j = 0; j < (int)gSize && head != NULL; d_j++, K_j++, j++){
		if (j < head->col){
			sum += (-(*K_j * k_i_M) * (*d_j));
		}
		else {
			sum += ((1 - (*K_j * k_i_M)) * (*d_j));
			head = head->next;
		}
	}

	while (j < (int)gSize)
	{
		sum += (-(*K_j * k_i_M) * (*d_j));
		j++;
		d_j++;
		K_j++;
	}

	dQ = (4 * d_i * sum) + (4 * k_i_M * k_i);
	/* Restore s to initial state */
	s[moved_v] *= -1;

	return dQ;
}



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
void add_after(Node *head, int i, double orig_val){
	head->next = (Node*)malloc(sizeof(Node));
	VERIFY(head->next != NULL,MEM_ALLOC_ERROR)
	head->next->col = i;
	head->next->val = orig_val;
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
			add_after(subHead, i, AHead->val);
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
	subHead->next = NULL; 
	return rowSize;
}

/* creating a sub sparse matrix for Algorithm 2 */
spmat* create_sub_sparse_matrix_linked(const struct _spmat  *A, Subgroup g, int sizeG){
	spmat *sub = spmat_allocate_list(sizeG);
	int_vector spmatSize = sub->spmatSize;
	Node **subMat = sub->private, **AMat = A->private, *AHead, *subHead, *tmp;
	Subgroup row;
	num k;
	
	for (row = g ; row < g + sizeG ; row++){
		subHead = (Node*)malloc(sizeof(Node));
		VERIFY(subHead != NULL,MEM_ALLOC_ERROR)
		tmp = subHead;
		AHead = AMat[*row];
		if (AHead != NULL){
			k = add_row_to_sub_linked(g, *row, g, AHead, subHead, sizeG);
		}
		else{
			k = 0;
			subHead->next = NULL;
		}
		*spmatSize++=k;
		subHead = subHead->next;
		free(tmp);
		*(subMat++) = subHead;
	}
	
	return sub;
}


void add_row_linked(spmat *A, const double *row, int i){
	Node *head = NULL, *tail;
	int j = 0, n = A->n;
	const DATA *rowPtr = row;
	while (j < n && !NON_ZERO(*rowPtr)){ /* getting to the first non-zero value and adding it to the first node*/
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
			while (j < n && !NON_ZERO(*rowPtr)){
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
	free(A->spmatSize);
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
	A->get_modularity_score = get_modularity_score_linked;
	A->create_sub_mat=create_sub_sparse_matrix_linked;
	A->spmatSize=(int_vector)malloc(n*sizeof(num));
	VERIFY(A->spmatSize!=NULL, MEM_ALLOC_ERROR)
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
	num *rp=imp->rowptr+i;
	scalar *valOffset=imp->values + *rp;
	num *colOffset=imp->colind + *rp;
	num rowNNZ=0, n=A->n;
	const double *row_p=row;

	for (row_p=row; row_p < row+n; row_p++)
		if (NON_ZERO(*row_p)){
			rowNNZ++;
			*valOffset++=*row_p;
			*colOffset++=row_p-row;
		}

	*(rp+1)=*rp+rowNNZ;
}

void get_row_arrays(const struct _spmat *A, int i, vector row){
	arraymat *imp=((arraymat*)(A->private));
	num j = *(imp->rowptr+i);
	num n=A->n, row_nnz = *(imp->rowptr+i+1)-j;
	num *src_col = imp->colind+j;
	double *src_val = imp->values+j;
	memset(row, 0, n*sizeof(double));
	while (row_nnz-->0)
		*(row+*src_col++)=*src_val++;
}


void free_arrays(struct _spmat *A)
{
	arraymat *imp=((arraymat*)(A->private));
	free(imp->values);
	free(imp->colind);
	free(imp->rowptr);
	free(A->spmatSize);
	free(A->private);
	free(A);
}

void mult_arrays(const struct _spmat *A, const double *v, double *result)
{
	arraymat *imp=((arraymat*)(A->private));
	num *rows = imp->rowptr, *r;
	num *j, rowNNZ, n=(num)A->n;
	scalar sum, *u;
	for (r = rows; r < rows+n; r++)
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

num get_g_row_nnz_arrays(arraymat *orig, int i, Subgroup g, int sizeG){
	num *r_ptr = orig->rowptr+i;
	num  r_i=*r_ptr, row_nnz=*(r_ptr+1)-r_i;
	num *cols=orig->colind, *orig_c=cols+r_i;
	num cnt=0;
	Subgroup p=g;
	while (orig_c < cols + r_i + row_nnz && p < g+sizeG){
		if 	(*p < *orig_c)
			p++;
		else if (*p > *orig_c)
			orig_c++;
		else {
			cnt++;
			p++;
			orig_c++;
		}
	}
	return cnt;
}

num get_g_nnz_arrays(const struct _spmat *A, Subgroup g, int sizeG){
	arraymat *orig=((arraymat*)(A->private));
	Subgroup p;
	num cnt=0;
	for (p=g; p<g+sizeG; p++)
		cnt+=get_g_row_nnz_arrays(orig, *p, g, sizeG);
	return cnt;
} 

num add_row_to_sub_arrays(arraymat *orig, int orig_i, arraymat *dest, int sub_i, Subgroup g, num sizeG){
	num 	k=orig->rowptr[orig_i], orig_nnz=orig->rowptr[orig_i+1]-k;
	num 	*c_orig=orig->colind+k;
	num 	*r_dest=dest->rowptr+sub_i;
	num		*c_dest=dest->colind, *c_dest_p;
	double 	*v_orig=orig->values+k;
	double 	*v_dest=dest->values;
	Subgroup p=g;
	
	if (sub_i==0)
		*r_dest=0;
	else {
		c_dest += *r_dest;
		v_dest += *r_dest;
	}

	c_dest_p=c_dest;
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
			*(c_dest_p++)=p-g;
			c_orig++;
			p++;
		}
	}

	*(r_dest+1)= *r_dest + (c_dest_p-c_dest); 
	return (c_dest_p-c_dest); 
}


spmat* create_sub_sparse_matrix_array(const struct _spmat  *A, Subgroup g, int sizeG){
	arraymat *orig=(arraymat*)(A->private), *imp;
	spmat *sub;
	Subgroup p;
	num sub_nnz = get_g_nnz_arrays(A, g, sizeG);
	
	sub=spmat_allocate_array(sizeG, sub_nnz);
	imp=(arraymat*)(sub->private);
	for (p=g; p<g+sizeG; p++)
		sub->spmatSize[p-g] = add_row_to_sub_arrays(orig, *p, imp, p-g, g, sizeG);
	
	return sub;
}

double get_modularity_score_array(const struct _spmat *A, vector s, int moved_v, int_vector K, num M){
	arraymat *mat = (arraymat*)(A->private);
	num gSize = A->n;
	double dQ;
	vector d_j;
	int_vector K_j;
	num k_i = K[moved_v];
	num j = mat->rowptr[moved_v], A_ij;
	num row_nnz = mat->rowptr[moved_v+1]-j;
	int_vector colind = mat->colind+j;
	int d_i;
	double k_i_M = k_i / (double)M, sum = 0;
	
	s[moved_v] *= -1;
	d_i = s[moved_v];

	
	for (d_j = s, K_j = K; K_j < K + gSize; d_j++, K_j++){
		A_ij = (K_j-K == *colind ? 1 : 0);
		row_nnz -= A_ij;
		sum += ((A_ij-(*K_j * k_i_M)) * (*d_j));
		colind += A_ij && row_nnz>0 ? 1 : 0;
	}

	dQ = (4 * d_i * sum) + (4 * k_i_M * k_i);
	/* Restore s to initial state */
	s[moved_v] *= -1;

	return dQ;
}


/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
spmat* spmat_allocate_array(int n, int nnz)
{
	spmat *A;
	double *values;
	int_vector colind, rowptr;

	A=(spmat*)malloc(sizeof(spmat));
	VERIFY(A != NULL, MEM_ALLOC_ERROR)

	A->spmatSize=(int_vector)malloc(n*sizeof(num));
	VERIFY(A->spmatSize!=NULL, MEM_ALLOC_ERROR)

	A->private=(arraymat*)malloc(sizeof(arraymat));
	VERIFY(A->private != NULL, MEM_ALLOC_ERROR)

	values = (scalar*)malloc(nnz * sizeof(scalar));
	VERIFY(values != NULL, MEM_ALLOC_ERROR)
	colind = (int_vector)malloc(nnz * sizeof(num));
	VERIFY(colind != NULL, MEM_ALLOC_ERROR)
	rowptr = (int_vector)malloc((n+1) * sizeof(num));
	VERIFY(rowptr != NULL, MEM_ALLOC_ERROR)
	*rowptr = 0; /*First non-zero element of matrix is always indexed as values[0], colind[0] */

	((arraymat*)A->private)->values = values;
	((arraymat*)A->private)->colind = colind;
	((arraymat*)A->private)->rowptr = rowptr;

	A->n=n;
	A->get_row=get_row_arrays;
	A->add_row=add_row_arrays;
	A->free=free_arrays;
	A->mult=mult_arrays;
	A->get_modularity_score = get_modularity_score_array;
	A->create_sub_mat=create_sub_sparse_matrix_array;
	return A;
}