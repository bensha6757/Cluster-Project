#ifndef LEADING_EIGENPAIR_H_
#define LEADING_EIGENPAIR_H_


#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Types.h"

#define IS_POSITIVE(X) ((X) > 0.00001)

/**
 * Compute leading eigen pair of modularity Matrix B_hat[g]
 */
void leading_eigenpair(modMat *B, modMat *Bg, vector leadEigenVec, double* leadEigenVal);

/**
 * Compute dot product of two vectors of length d
 */
double dot_prod(vector v, vector u, size_t d);

#endif /* LEADING_EIGENPAIR_H_ */