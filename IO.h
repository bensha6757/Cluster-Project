#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include "Stack.h"
#include "io_mem_errors.h"
#include "modMat.h"
#include "spmat.h"
#include "Types.h"

/** Loads an input file into a modmat struct in memory.
 * @param filename - a string, path of the input file.
 * @param mat - an address for the modmat struct to be stored at. 
 * */ 
void load_input_file(char* filename, modMat **mat);

void generate_output_file(Stack *O, char *outputPath);

#endif /* IO_H_ */