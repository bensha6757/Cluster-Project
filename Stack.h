typedef int bool;

typedef struct Snode{
	int* g;
	int sizeG;
	Snode* next;
} Snode;

typedef struct Stack{
	int size;
	Snode* top;
} Stack;

void init(Stack* s);

void push(Stack* s, int* g, int sizeG);

int* pop(Stack* s, int* sizeG);

bool isEmpty(Stack* s);

void delete_Stack(Stack* s);

void delete_StackNodes(Snode* node);