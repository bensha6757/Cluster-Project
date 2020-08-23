#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include "Assertions.h"
#include "Stack.h"
#include "modMat.c"

void add_to_stacks(Stack* P, Stack* O, int sizeG, int sizeG1, int sizeG2, int* g1, int* g2, int* g){
	if (sizeG1 == 0 || sizeG2 == 0){
			push(O, g, sizeG);
	}
	else{
		free(g); /* g has been divided into 2 groups, so now it has no use anymore */
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
/* Algorithm 3*/
Stack* div_into_mod_groups(modMat* B, int *g, int sizeG){
	Stack* P = (Stack*)malloc(sizeof(Stack)), *O = (Stack*)malloc(sizeof(Stack));
	Snode* node;
	int *res, sizeG1, sizeG2, *g1, *g2;
	verify(P != NULL);
	verify(O != NULL);
	init(P);
	init(O); /* initializing O to be empty */
	push(P,g,sizeG); /* starting with a trivial division into one group */
	while (!isEmpty(P)){
		g = pop(P, &sizeG); /* g would be g1 or g2 from some previous iteration */
		div_into_two(B, g, sizeG, &g1, &g2, &sizeG1, &sizeG2); /* calling Algorithm 2, (g1,g2,sizeG1,sizeG2) will be populate inside the function with the result of the division */
		add_to_stacks(P, O, sizeG, sizeG1, sizeG2, g1, g2, g);
	}
	free(P);
	return O;
}