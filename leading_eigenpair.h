#ifndef LEADING_EIGENPAIR_H_
#define LEADING_EIGENPAIR_H_


#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Types.h"

#define EPSILON 0.00001
#define IS_POSITIVE(X) ((X) > (EPSILON))
#define LIMIT_ITER(N) ((0.5)*(N)*(N) + 10000*(N)+30000)


/**
 * Compute dot product of two vectors of length d.
 */
double dot_prod(vector v, vector u, num d);

/** Compute leading eigen pair of modularity Matrix B_hat[g], and store results in referenced pointers.
 * @param Bg - a modularity matrix.
 * @param leadEigenVec - the address in memory to store the leading eigenvector (vector) in.
 * @param leadEigenVal - the address in memory to store the leading eigenvalue (scalar) in.
 */
scalar leading_eigenpair(modMat *Bg, vector *leadEigenVec);


#endif /* LEADING_EIGENPAIR_H_ */