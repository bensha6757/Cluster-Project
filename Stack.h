typedef int bool;
typedef size_t* Subgroup;

typedef struct Snode{
	Subgroup g;
	size_t sizeG;
	Snode* next;
} Snode;

typedef struct Stack{
	size_t size;
	Snode* top;
} Stack;

void init(Stack* s);

void push(Stack* s, Subgroup g, size_t sizeG);

Subgroup pop(Stack* s, size_t* sizeG);

bool isEmpty(Stack* s);

void delete_Stack(Stack* s);