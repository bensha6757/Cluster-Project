#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"
#include "Stack.h"
#include "Types.h"

void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat=NULL;
	Subgroup g, p;
	Stack *O;
	load_input_file(inputFileName,&mat); /* reading input */
	g = (Subgroup)malloc(mat->gSize * sizeof(num));
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	for (p=g; p<g+mat->gSize; p++)
		*p=p-g;
	O = div_into_mod_groups(mat, g, mat->gSize); /* calling Algorithm 3 */
	mat->free(mat);
	/*free(g);*/
	generate_output_file(O, outputFileName); /* writing result to file */
	delete_Stack(O);
	#ifdef DEBUG
	printf("SUCCESS: runClusterProject\n");
	#endif
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	clock_t start, end;

	#ifdef DEBUG
	printf("Received file names: %s, %s\n",inputFileName,outputFileName);
	#endif
	VERIFY(argc - 1 == 2, MISSING_ARG_ERROR)
	srand(time(NULL));
	start = clock();
	runClusterProject(inputFileName, outputFileName);
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	return SUCCESS;
}