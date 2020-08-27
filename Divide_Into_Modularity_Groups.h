#ifndef DIVIDE_INTO_MODULARITY_GROUPS_H_
#define DIVIDE_INTO_MODULARITY_GROUPS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include "io_mem_errors.h"
#include "Types.h"
#include "Stack.h"
#include "modMat.h"
#include "divide_into_two.h"


Stack* div_into_mod_groups(modMat* B, Subgroup g, size_t sizeG);


#endif /* DIVIDE_INTO_MODULARITY_GROUPS_H_ */