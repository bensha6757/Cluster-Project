#include "Divide_Into_Modularity_Groups.h"
#include "IO.h"
#include "Stack.h"
#include "Types.h"

void runClusterProject(char* inputFileName, char* outputFileName){
	modMat *mat=NULL;
	Subgroup g, g_i;
	Stack *O;
	num i, gSize;
	load_input_file(inputFileName,&mat); /* reading input */
	
	gSize = mat->gSize;	
	g = (Subgroup)malloc(gSize * sizeof(num)); /* generating trivial division into one group */
	VERIFY(g!=NULL,MEM_ALLOC_ERROR)
	for (g_i = g , i = 0 ; i < gSize ; g_i++, i++){
		*g_i = i;
	}
	
	O = divide_into_mod_groups(mat, g, gSize); /* calling Algorithm 3 */

	mat->free(mat); 
	/*free(g);*/

	generate_output_file(O, outputFileName); /* writing result to file */

	delete_Stack(O);
	
}

int main(int argc, char* argv[]){
	char* inputFileName=argv[1];
	char* outputFileName=argv[2];
	clock_t start, end;

	VERIFY(argc - 1 == 2, MISSING_ARG_ERROR)
	srand(time(NULL));
	start = clock();
	runClusterProject(inputFileName, outputFileName);
	end = clock();
	printf("Execution took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	return SUCCESS;
}