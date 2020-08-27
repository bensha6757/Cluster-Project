#ifndef STACK_H_
#define STACK_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "modMat.h"
#include "Types.h"

typedef struct Snode{
	Subgroup g;
	size_t sizeG;
	struct Snode* next;
} Snode;

typedef struct Stack{
	size_t size;
	Snode* top;
} Stack;

void init(Stack* s);

void push(Stack* s, Subgroup g, size_t sizeG);

Subgroup pop(Stack* s, size_t* sizeG);

boolean isEmpty(Stack* s);

void delete_Stack(Stack* s);

#endif /* STACK_H_ */