#include <stdio.h>
#include "Stack.h"
#include "io_mem_errors.h"
#include "modMat.h"
#include "spmat.h"

void generate_input_file(char* filename, modMat *mat);

void generate_output_file(Stack *O, char *outputPath);

void convert_adj_list(size_t k, size_t n, int_vector adj, double *res);

size_t read_totalV_from_file(FILE *input);

void load_mod_matrix_from_file(FILE *input, modMat *B);