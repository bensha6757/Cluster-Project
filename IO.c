#include "IO.h"

/*********************************
 * INPUT FILE TO MEMORY FUNCTIONS*
 *********************************/


/*	Converts node indices (like input file format) to a {0,1} vector of length n
 * 	that can be added to a spmat via add_row. 
 */

vector convert_adj_list(int_vector adj, num adj_len, num n){
	int_vector k;
	vector line;
	line=(double*)calloc(n,sizeof(double));
	VERIFY(line!=NULL, MEM_ALLOC_ERROR)
	for (k=adj; k<adj+adj_len; k++)
		*(line+*k)=1;
	return line;
}



num read_network_size_from_file(FILE *input, num *edges_num){
	num n, i, k, *tmp, cnt=0;
	VERIFY(fread(&n,sizeof(int),1,input) == 1,FILE_READ_ERROR)
	tmp=(int_vector)malloc(n*sizeof(num));
	VERIFY(tmp!=NULL,MEM_ALLOC_ERROR)
	/* Count number of edges (for spmat arrays impl.) */
	for (i=0; i<n; i++){
		VERIFY(fread(&k,sizeof(int),1,input) == 1,FILE_READ_ERROR)
		VERIFY(fread(tmp,sizeof(int),k,input) == k,FILE_READ_ERROR)
		cnt+=k;
	}
	*edges_num=cnt/2; /*Assuming sum_v(deg(v))==2*|E| */
	free(tmp);
	rewind(input);
	return n;
}

void load_mod_matrix_from_file(FILE *input, modMat *B){
	int_vector K=B->K, g=B->g, sp=B->spmatSize, neighbours;
	num k_i, i;
	vector matLine=NULL;
	VERIFY(fread(&i,sizeof(num),1,input) == 1, FILE_READ_ERROR)  /*Assuming file rewinded, skip |V| */
	/* Populate B with A, K matrices, and compute M */
	for (i=0; i<B->gSize; i++){
		VERIFY(fread(&k_i,sizeof(num),1,input) == 1, FILE_READ_ERROR)
		neighbours=(int_vector)malloc(k_i*sizeof(num));
		VERIFY(neighbours!=NULL, MEM_ALLOC_ERROR)
		VERIFY(fread(neighbours,sizeof(num), k_i, input) == k_i, FILE_READ_ERROR)
		matLine=convert_adj_list(neighbours, k_i, B->gSize);
		B->A->add_row(B->A,matLine,i);
		*(g++)=i; /*Fill g subgroup array with 0,1,2,...,n */
		*(K++)=k_i;
		*(sp++)=k_i;
		B->M+=k_i;
		free(neighbours);
		free(matLine);
	}
	B->currM=B->M;
	rewind(input);
	set_1_norm(B);
}

/*** INTERFACE FOR MAIN PROGRAM ***/

void load_input_file(char* filename, modMat **mat){
	num N, M;
	FILE* inputFile = fopen(filename,"r");
	N=read_network_size_from_file(inputFile, &M);
	*mat = allocate_mod_mat(N);
	load_mod_matrix_from_file(inputFile, *mat);
	fclose(inputFile);
}

void generate_output_file(Stack *O, char *outputPath){
    FILE *out = fopen(outputPath, "w"); /*file should be open for write in main*/
    Snode *node = O->top;
    VERIFY(out != NULL, FILE_WRITE_ERROR)
	VERIFY(fwrite(&(O->size), sizeof(num), 1, out) == 1, FILE_WRITE_ERROR)
    while(node != NULL){
        VERIFY(fwrite(&(node->sizeG),sizeof(num),1, out) == 1, FILE_WRITE_ERROR)
        VERIFY(fwrite(node->g, sizeof(num), node->sizeG, out) == node->sizeG , FILE_WRITE_ERROR)
        node = node->next;
    }
    fclose(out);
}