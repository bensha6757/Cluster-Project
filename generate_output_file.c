#include <stdio.h>
#include "Stack.h"
#include "Assertions.h"


void generate_output_file(Stack *O, char *outputPath){
    FILE *out = fopen(outputPath, "w"); /*file should be open for write in main*/
    Snode *node = O->top;
    verify(out != NULL);
    fputc(O->size + '0', out);
    while(node != NULL){
        fputc(node->sizeG + '0', out);
        for (int *gi = node->g ; gi - node->g < node->sizeG ; gi++){
            fputc(*gi + '0', out);
        }
        node = node->next;
    }
    fclose(out);
}
