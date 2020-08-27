#include "Stack.h"

void init(Stack *s){
	s->size = 0;
	s->top = NULL;
}

/* adding the new node to the beggining of the stack (implemented as linked-list), and updating the stack size*/
void push(Stack *s, Subgroup g, size_t sizeG){ 
	Snode * node = (Snode*)malloc(sizeof(Snode));
	VERIFY(node != NULL,MEM_ALLOC_ERROR)
	node->g = g;
	node->sizeG = sizeG;
	node->next = s->top;
	s->top = node;
	s->size++;
}

Subgroup pop(Stack *s, size_t *sizeG){
	Subgroup g = s->top->g;
	Snode * p = s->top;
	*sizeG = p->sizeG;
	s->top = s->top->next;
	s->size--;
	free(p);
	return g;
}

boolean isEmpty(Stack *s){
	return s->size == 0;
} 

void delete_StackNodes(Snode *node){
	if (node != NULL){
		delete_StackNodes(node->next);
		free(node->g);
		free(node);
	}
}

void delete_Stack(Stack *s){
	delete_StackNodes(s->top);
	free(s);
}

