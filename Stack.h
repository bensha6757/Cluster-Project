#ifndef STACK_H_
#define STACK_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "modMat.h"
#include "Types.h"

/****************************************************************************************************************
 * A data structure storing subgroups (communities) of vertices for the modularity division algorithm (Alg. 3).
 * Allows a single Subgroup to be pushed or popped at request. 
 ****************************************************************************************************************/

typedef struct Snode{
	Subgroup g;
	num sizeG;
	struct Snode* next;
} Snode;

typedef struct Stack{
	num size;
	Snode* top;
} Stack;

void init(Stack* s);

/* adding the new node to the beggining of the stack (implemented as linked-list),
** and updating the stack size
*/
void push(Stack* s, Subgroup g, num sizeG);

/* pop the first element from the stack, and free the popped element.
** returning the subgroup that was popped 
*/
Subgroup pop(Stack* s, num *sizeG);

boolean isEmpty(Stack* s);

void delete_Stack(Stack* s);

#endif /* STACK_H_ */