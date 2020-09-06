#ifndef LINKED_LIST_H_
#define LINKED_LIST_H_

#include "io_mem_errors.h"
#include "Types.h"
#include <stdlib.h>

typedef struct mNode{
	num ind;
	struct mNode* next;
} mNode;

typedef struct Linked_list_moved{
	mNode *head;
	num size;
} Linked_list_moved;


Linked_list_moved *create_unmoved(num size);

void delete_node(Linked_list_moved *l, mNode *prev, mNode *curr);

void delete_unmoved(Linked_list_moved *s);

#endif /* LINKED_LIST_H_ */