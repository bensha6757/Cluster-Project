#ifndef DIVIDE_INTO_MODULARITY_GROUPS_H_
#define DIVIDE_INTO_MODULARITY_GROUPS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include "IO_Mem_Errors.h"
#include "Types.h"
#include "Stack.h"
#include "modMat.h"
#include "Divide_Into_Two.h"

/** Divides a network represented as modularity matrix B and a subgroup of indices g (with number of indices sizeG) to
 *  communities. Implements Algorithm 3 of the project.
 *  Returns a stack of communinties, represented as sorted Subgroup object of indices.
 */
Stack* divide_into_mod_groups(modMat* B, Subgroup g, num sizeG);


#endif /* DIVIDE_INTO_MODULARITY_GROUPS_H_ */