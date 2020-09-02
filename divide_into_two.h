#ifndef DIVIDE_INTO_TWO_H_
#define DIVIDE_INTO_TWO_H_

#include "io_mem_errors.h"
#include "modMat.h"
#include "leading_eigenpair.h"
#include "Types.h"
#include <stdio.h>
#include <stdlib.h>
typedef enum div_res_t {
	GROUP_INDIVISIBLE,
	GROUP_DIVIDED
} DIV_RESULT;

DIV_RESULT div_into_two(modMat *B,Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2);

#endif /* DIV_INTO_TWO_H_ */