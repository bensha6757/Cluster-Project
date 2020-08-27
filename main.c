#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"
#include "Stack.h"
#include "Types.h"


void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat=NULL;
	Subgroup g;
	Stack *O;
	load_input_file(inputFileName,mat); /* reading input */
	g = (Subgroup)malloc(mat->gSize * sizeof(size_t));
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	memcpy(g, mat->g, mat->gSize * sizeof(size_t));
	O = div_into_mod_groups(mat, g, mat->gSize); /* calling Algorithm 3 */
	generate_output_file(O, outputFileName); /* writing result to file */
	mat->free(mat);
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	VERIFY(argc == 2,FILE_READ_ERROR)
	runClusterProject(inputFileName, outputFileName);
	return 0;
}