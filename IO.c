#include "IO.h"
#define DEBUG

/*********************************
 * INPUT FILE TO MEMORY FUNCTIONS*
 *********************************/


/*	Convert node indices (like input file format) to a pre-allocated result {0,1} vector of length n
 * 	that can be added to a spmat via add_row. 
 */
void convert_adj_list(size_t deg, size_t n, int* adj, double *res){
	int *i, temp;
	for (i=adj; i<adj+deg; i++){
		temp=*i;
		while (temp-->0)
			*res++=0;
		*res++=1;
	}
	for(temp=n-*i; temp>0; temp--)
		*res++=0;
}


size_t read_network_size_from_file(FILE *input){
	int n;
	VERIFY(fread(&n,sizeof(size_t),1,input) == 1,FILE_READ_ERROR)
	rewind(input);
	return (size_t)n;
}

void load_mod_matrix_from_file(FILE *input, modMat *B){
	int i;
	int_vector K=B->K, g=B->g;
	int *neighbours, k_i;
	double *matLine=(double*)malloc(B->gSize*sizeof(double));
	VERIFY(matLine!=NULL, MEM_ALLOC_ERROR)
	VERIFY(fread(&i,sizeof(size_t),1,input) == 1, FILE_READ_ERROR)  /*Assuming file rewinded, skip |V| */
	/* Populate B with A, K matrices, and compute M */
	for (i=0; i<(int)B->gSize; i++){
		VERIFY(fread(&k_i,sizeof(int),1,input) == 1, FILE_READ_ERROR)
		*(K++)=(size_t)k_i;
		B->M+=(size_t)k_i;
		neighbours=(int*)malloc(k_i*sizeof(int));
		VERIFY(neighbours!=NULL, MEM_ALLOC_ERROR)
		VERIFY(fread(neighbours,sizeof(size_t),k_i,input) == (size_t)k_i, FILE_READ_ERROR)
		convert_adj_list((size_t)k_i, B->gSize, neighbours, matLine);
		B->A->add_row(B->A,matLine,i);
		*(g++)=i; /*Fill g subgroup array with 0,1,2,...n */
		free(neighbours);
	}
	free(matLine);
	rewind(input);
	set_1_norm(B);
}

/*** INTERFACE FOR MAIN PROGRAM FUNCTIONS ***/

void load_input_file(char* filename, modMat *mat){
	size_t N;
	FILE* inputFile = fopen(filename,"r");
	#ifdef DEBUG
	printf("Open input file\n");
	#endif
	N=read_network_size_from_file(inputFile);
	#ifdef DEBUG
	printf("Open Read |V|=%d\n", (int)N);
	#endif
	mat = allocate_mod_mat(N);
	load_mod_matrix_from_file(inputFile, mat);
	fclose(inputFile);
}

void generate_output_file(Stack *O, char *outputPath){
    FILE *out = fopen(outputPath, "w"); /*file should be open for write in main*/
    Snode *node = O->top;
	Subgroup gi;
    VERIFY(out != NULL, FILE_WRITE_ERROR)
    fputc(O->size + '0', out);
    while(node != NULL){
        fputc(node->sizeG + '0', out);
        for (gi = node->g ; gi  < node->sizeG + node->g ; gi++){
            fputc(*gi + '0', out);
        }
        node = node->next;
    }
    fclose(out);
}