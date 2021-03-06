#include "IO.h"

/*********************************
 * INPUT FILE TO MEMORY FUNCTIONS*
 *********************************/


/**	Converts adj array of length adj_len, consisting of node indices (like in input file format) 
 *  to a preallocated {0,1} vector of length n that can be added to an instance of spmat via add_row. 
 */
void convert_adj_list(int_vector adj, num adj_len, vector line, num n){
	int_vector k;
	memset(line, 0, n*sizeof(double));
	for (k=adj; k<adj+adj_len; k++)
		*(line+*k)=1;
}

/** Returns the number of nodes in network (i.e. n == |V|) represented as a formatted file.
 * 	Additionaly, computes sum the of degrees of the network (i.e. 2*|E|).
 *  @param input - a previosly open input file, formatted as a network's adjecancy list.
 * 	@param m - an address of non-negative integer to store sum of network's degrees. 
 */
num read_network_size_and_nnz_from_file(FILE *input, num *m){
	num n, i, k;
	num* tmp;
	VERIFY(fread(&n,sizeof(num),1,input) == 1,FILE_READ_ERROR)
	tmp=(num*)calloc(n,sizeof(num));
	VERIFY(tmp!=NULL, MEM_ALLOC_ERROR)
	for (i=0; i<n; i++){
		VERIFY(fread(&k,sizeof(num),1,input) == 1,FILE_READ_ERROR)
		*m+=k;
		VERIFY(fread(tmp,sizeof(num),k,input) == k,FILE_READ_ERROR)
	}
	free(tmp);
	rewind(input);
	return n;
}

/** Reads an input file into a modularity matrix in memory.
 * 	Assumes file will no longer be used afterwards.
 *  @param input - a previosly open input file, formatted as a network's adjecancy list.
 * 	@param B - a previously allocated modularity matrix B object.
 */
void load_mod_matrix_from_file(FILE *input, modMat *B){
	int_vector K = B->K, sp=B->A->spmatSize, neighbours;
	num k_i, i, gSize = B->gSize, M = 0;
	vector matLine;
	matLine = (double*)malloc(gSize * sizeof(double));
	VERIFY(matLine!=NULL, MEM_ALLOC_ERROR);
	neighbours = (int_vector)malloc(gSize * sizeof(num));
	VERIFY(neighbours!=NULL, MEM_ALLOC_ERROR)
	VERIFY(fread(&i,sizeof(num),1,input) == 1, FILE_READ_ERROR) /*Assume reading the file after rewind */
	/* Populate B with A, K matrices, and compute M */
	for (i = 0 ; i < gSize ; i++){
		VERIFY(fread(&k_i,sizeof(num),1,input) == 1, FILE_READ_ERROR)
		VERIFY(fread(neighbours, sizeof(num), k_i, input) == k_i, FILE_READ_ERROR)
		convert_adj_list(neighbours, k_i, matLine, gSize);
		B->A->add_row(B->A, matLine, i);
		*(K++) = k_i;
		*(sp++) = k_i;
		M += k_i;
	}
	free(neighbours);
	free(matLine);
	B->M = M;
	B->currM = M;
	set_1_norm(B);
}

/*** INTERFACE FOR MAIN PROGRAM ***/


void load_input_file(char* filename, modMat **mat){
	num N, M=0;
	FILE* inputFile;
	inputFile = fopen(filename,"r");
	VERIFY(inputFile != NULL, FILE_NOT_FOUND_ERROR)
	N = read_network_size_and_nnz_from_file(inputFile, &M);
	*mat = allocate_mod_mat(N, M, FALSE);
	VERIFY(*mat != NULL, MEM_ALLOC_ERROR)
	load_mod_matrix_from_file(inputFile, *mat);
	fclose(inputFile);
}

void generate_output_file(Stack *O, char *outputPath){
    FILE *out = fopen(outputPath, "w"); /*file should be open for write in main*/
    Snode *node = O->top;
    VERIFY(out != NULL, FILE_NOT_FOUND_ERROR)
	VERIFY(fwrite(&(O->size), sizeof(num), 1, out) == 1, FILE_WRITE_ERROR)
    while(node != NULL){
        VERIFY(fwrite(&(node->sizeG),sizeof(num),1, out) == 1, FILE_WRITE_ERROR)
        VERIFY(fwrite(node->g, sizeof(num), node->sizeG, out) == node->sizeG , FILE_WRITE_ERROR)
        node = node->next;
    }
    fclose(out);
}