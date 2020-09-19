#include <stdio.h>
#include <stdlib.h>

#define VERIFY(expression)                                 \
    if (!(expression)){                                     \
        printf("Verification Failed.\n");                    \
        exit(1);                                              \
    }         

void printList(unsigned int *list, unsigned int len){
    unsigned int *i;
    for (i=list; i<list+len; i++)
        printf("%d, ",(int)(*i));
}

void printBinaryFile(FILE *f){
    unsigned int i, n, k, *list;
    VERIFY(fread(&n,sizeof(unsigned int),1,f)==1)
    printf("# Groups: %d\n",(int)n);
    for (i=0; i<n; i++){
        VERIFY(fread(&k,sizeof(unsigned int),1,f)==1)
        printf("Group %d of size %d: ", i+1, k);
        list=(unsigned int*)malloc(k*sizeof(unsigned int));
        VERIFY(list!=NULL);
        VERIFY(fread(list,sizeof(unsigned int),k,f)==k)
        printList(list,k);
        printf("\n");
    }
}

int main(int argc, char* argv[]){
	FILE *f;
	VERIFY(argc - 1 == 1)
    f=fopen(argv[1],"r");
	printBinaryFile(f);
    fclose(f);
	return 0;
}