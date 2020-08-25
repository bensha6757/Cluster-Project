#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"


void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat;
	Subgroup g;
	Stack *O;
	load_input_file(inputFileName,mat);
	g = (Subgroup)malloc(mat->gSize * sizeof(size_t));
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	memcpy(g, mat->g, mat->gSize * sizeof(size_t));
	O = div_into_mod_groups(mat, g, mat->gSize);
	generate_output_file(O, outputFileName);
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	runClusterProject(inputFileName, outputFileName);
	return 0;
}