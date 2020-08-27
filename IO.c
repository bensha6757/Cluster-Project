#include "IO.h"


/*********************************
 * INPUT FILE TO MEMORY FUNCTIONS*
 *********************************/

/*  Read bytes from file or exit if failed. used to replace assert functionality */
void read_num_from_file(int *dest,unsigned int num,unsigned int size,FILE* file){
	VERIFY(fread(dest,num,size,file) == size,FILE_READ_ERROR)
}

/*	Convert node indices (like input file format) to a pre-allocated result {0,1} vector of length n
 * 	that can be added to a spmat via add_row. */
void convert_adj_list(size_t k, size_t n, int* adj, double *res){
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


size_t read_totalV_from_file(FILE *input){
	int n;
	read_num_from_file(&n,sizeof(size_t),1,input);
	rewind(input);
	return n;
}

void load_mod_matrix_from_file(FILE *input, modMat *B){
	int i;
	int_vector k=B->K, g=B->g;
	int *inputNeighbours, currDeg;
	double *matLine=(double*)malloc(B->gSize*sizeof(size_t));
	VERIFY(matLine!=NULL,MEM_ALLOC_ERROR)
	read_num_from_file(&i,sizeof(size_t),1,input); /*Assuming file rewinded, skip |V| */
	/* Populate B with A, K matrices, and compute M */
	for (i=0; i<(int)B->gSize; i++){
		read_num_from_file(&currDeg,sizeof(int),1,input);
		*(k++)=(size_t)currDeg;
		B->M+=(size_t)currDeg;
		inputNeighbours=(int*)malloc(currDeg*sizeof(int));
		VERIFY(inputNeighbours!=NULL,MEM_ALLOC_ERROR)
		read_num_from_file(inputNeighbours,sizeof(size_t),currDeg,input);
		convert_adj_list((size_t)currDeg, B->gSize, inputNeighbours, matLine);
		B->A->add_row(B->A,matLine,i);
		*(g++)=i; /*Fill g subgroup array with 0,1,2,...n */
		free(inputNeighbours);
	}
	free(matLine);
	rewind(input);
	set_1_norm(B);
}

/*** INTERFACE FOR MAIN PROGRAM FUNCTIONS ***/

void load_input_file(char* filename, modMat *mat){
	size_t N;
	FILE* inputFile = fopen(filename,"r");
	N=read_totalV_from_file(inputFile);
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
        for (gi = node->g ; gi  < node->sizeG +  node->g ; gi++){
            fputc(*gi + '0', out);
        }
        node = node->next;
    }
    fclose(out);
}