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


Stack* divide_into_mod_groups(modMat* B, Subgroup g, num sizeG);


#endif /* DIVIDE_INTO_MODULARITY_GROUPS_H_ */