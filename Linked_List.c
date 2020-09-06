#include "Linked_List.h"

Linked_list_moved *create_unmoved(num size){
    num i;
    mNode *head, *tail;
    Linked_list_moved *l = (Linked_list_moved*)malloc(sizeof(Linked_list_moved));
    VERIFY(l!=NULL,MEM_ALLOC_ERROR)
    l->head = NULL;
    l->size = size;
    if (size == 0)
        return l;
    head = (mNode*)malloc(sizeof(mNode));
    VERIFY(head!=NULL,MEM_ALLOC_ERROR)
    head->ind = 0;
    tail = head;
    for (i = 1; i < size; i++){
        tail->next = (mNode*)malloc(sizeof(mNode));
        VERIFY(tail->next!=NULL,MEM_ALLOC_ERROR)
        tail->next->ind = i;
        tail = tail->next;
    }
    tail->next = NULL;
    l->head = head;
    return l;
}

void delete_node(Linked_list_moved *l, mNode *prev, mNode *curr){
    if(prev == NULL){
        l->head = curr;
        return;
    }
    prev->next = curr->next;
    free(curr);
}

void delete_list_nodes(mNode *node){
	if (node != NULL){
		delete_list_nodes(node->next);
		free(node);
	}
}

void delete_unmoved(Linked_list_moved *l){
	delete_list_nodes(l->head);
	free(l);
}