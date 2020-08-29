#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include "Stack.h"
#include "io_mem_errors.h"
#include "modMat.h"
#include "spmat.h"
#include "Types.h"


void load_input_file(char* filename, modMat *mat);

void generate_output_file(Stack *O, char *outputPath);

num read_totalV_from_file(FILE *input);

void load_mod_matrix_from_file(FILE *input, modMat *B);

#endif /* IO_H_ */