#ifndef IO_H_
#define IO_H_

#include <stdio.h>
#include "Stack.h"
#include "IO_Mem_Errors.h"
#include "modMat.h"
#include "Spmat.h"
#include "Types.h"

/**************************************************************************************
 * a module in charge of reading the input binary file and writing the final results, *
 *  containing the modularity division to groups                                      *
 * ************************************************************************************/

/** Loads an input file into a modularity matrix struct in memory.
 * @param filename - a string, path of the input file.
 * @param mat - an address for the modmat struct to be stored at. 
 * */ 
void load_input_file(char* filename, modMat **mat);

/** Unloads a stack of subgroups (communities) into a file.
 * @param O - A stack containing disjoint communities of a network.
 * @param outputPath - a string, path of the output file.
 * */ 
void generate_output_file(Stack *O, char *outputPath);

#endif /* IO_H_ */