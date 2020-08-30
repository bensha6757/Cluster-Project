
#include "Divide_Into_Modularity_Groups.h"

void add_to_stacks(Stack* P, Stack* O, num sizeG, num sizeG1, num sizeG2, Subgroup g1, Subgroup g2, Subgroup g){
	if (sizeG1 == 0 || sizeG2 == 0){
			push(O, g, sizeG);
	}
	else{
		if (sizeG1 == 1){
			push(O, g1, sizeG1);
		}
		else{
			push(P, g1, sizeG1);
		}
		if (sizeG2 == 1){
			push(O, g2, sizeG2);
		}
		else{
			push(P, g2, sizeG2);
		}
	}
}


void printG(Subgroup g, int n){
	Subgroup p;
	for (p = g; p < g + n ; p++){
		printf("%d, ", *p);
	}
	printf("\n");
}


/* Algorithm 3*/
Stack* div_into_mod_groups(modMat* B, Subgroup g, num sizeG){
	Stack* P = (Stack*)malloc(sizeof(Stack)), *O = (Stack*)malloc(sizeof(Stack));
	Subgroup g1, g2;
	num sizeG1, sizeG2;
	#ifdef DEBUG
	printf("BEGIN: div_into_mod_groups\n");
	#endif
	VERIFY(P != NULL,MEM_ALLOC_ERROR)
	VERIFY(O != NULL,MEM_ALLOC_ERROR)
	init(P);
	init(O); /* initializing O to be empty */
	push(P,g,sizeG); /* starting with a trivial division into one group */
	while (!isEmpty(P)){
		#ifdef DEBUG
		printf("Alg 3  - NEW ITER\n");
		#endif
		g = pop(P, &sizeG); /* g would be g1 or g2 from some previous iteration */
		div_into_two(B, g, sizeG, &g1, &g2, &sizeG1, &sizeG2); /* calling Algorithm 2, (g1,g2,sizeG1,sizeG2) will be populate inside the function with the result of the division */
		#ifdef DEBUG
		printf("Alg 2 completed - division:\n");
		printG(g1,sizeG1);
		printG(g2,sizeG2);
		#endif
		add_to_stacks(P, O, sizeG, sizeG1, sizeG2, g1, g2, g);
	}
	free(P);
	#ifdef DEBUG
	printf("SUCCESS: div_into_mod_groups\n");
	#endif
	return O;
}
