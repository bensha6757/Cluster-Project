#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"
#include "Stack.h"
#include "Types.h"

void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat=NULL;
	Subgroup g;
	Stack *O;
	load_input_file(inputFileName,&mat); /* reading input */
	#ifdef DEBUG
	printf("SUCCESS: load_input_file\n");
	#endif
	g = (Subgroup)malloc(mat->gSize * sizeof(num));
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	memcpy(g, mat->g, mat->gSize * sizeof(num));
	O = div_into_mod_groups(mat, g, mat->gSize); /* calling Algorithm 3 */
	#ifdef DEBUG
	printf("SUCCESS: div_into_mod g\n");
	#endif
	mat->free(mat);
	generate_output_file(O, outputFileName); /* writing result to file */
	delete_Stack(O);
	#ifdef DEBUG
	printf("SUCCESS: runClusterProject\n");
	#endif
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	#ifdef DEBUG
	printf("Received file names: %s, %s\n",inputFileName,outputFileName);
	#endif
	VERIFY(argc - 1 == 2, FILE_READ_ERROR)
	runClusterProject(inputFileName, outputFileName);
	return 0;
}