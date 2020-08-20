#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "Stack.h"

typedef struct Snode{
	int* g;
	int sizeG;
	Snode* next;
} Snode;

typedef struct Stack{
	int size;
	Snode* top;
} Stack;

void init(Stack *s){
	s->size = 0;
	s->top = NULL;
}

/* adding the new node to the beggining of the stack (implemented as linked-list), and updating the stack size*/
void push(Stack *s, int* g, int sizeG){ 
	Snode * node = (Snode*)malloc(sizeof(Snode));
	verify(node != NULL);
	node->g = g;
	node->sizeG = sizeG;
	node->next = s->top;
	s->top = node;
	s->size++;
}

int* pop(Stack *s, int *sizeG){
	int* g = s->top->g;
	Snode * p = s->top;
	*sizeG = p->sizeG;
	s->top = s->top->next;
	s->size--;
	free(p);
	return g;
}

bool isEmpty(Stack *s){
	return s->size == 0;
} 

void delete_Stack(Stack *s){
	delete_StackNodes(s->top);
	free(s);
}

void delete_StackNodes(Snode *node){
	if (node != NULL){
		delete_StackNodes(node->next);
		free(node->g);
		free(node);
	}
}