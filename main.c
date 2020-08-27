#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"
#include "Stack.h"
#include "Types.h"
#define DEBUG

void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat=NULL;
	Subgroup g;
	Stack *O;
	#ifdef DEBUG
	printf("Read input file begins");
	#endif
	load_input_file(inputFileName,mat); /* reading input */
	#ifdef DEBUG
	printf("Read input file successfully");
	#endif
	g = (Subgroup)malloc(mat->gSize * sizeof(size_t));
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	memcpy(g, mat->g, mat->gSize * sizeof(size_t));
	O = div_into_mod_groups(mat, g, mat->gSize); /* calling Algorithm 3 */
	mat->free(mat);
	generate_output_file(O, outputFileName); /* writing result to file */
	delete_Stack(O);
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	#ifdef DEBUG
	printf("Recieved file names: %s, %s",inputFileName,outputFileName);
	#endif
	VERIFY(argc - 1 == 2,FILE_READ_ERROR)
	runClusterProject(inputFileName, outputFileName);
	return 0;
}