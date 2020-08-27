#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include "io_mem_errors.h"
#include "Stack.h"
#include "modMat.h"


Stack* div_into_mod_groups(modMat* B, Subgroup g, size_t sizeG);