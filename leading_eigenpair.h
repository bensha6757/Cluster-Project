#ifndef LEADING_EIGENPAIR_H_
#define LEADING_EIGENPAIR_H_


#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Types.h"

 /*There are about O(N^2) Power iterations for a matrix of size N */
#define ITER_LIMIT(N) (0.5*(N)*(N)+10000*(N)+300000)


/**
 * Compute dot product of two vectors of length d.
 */
double dot_prod(vector v, vector u, num d);

/** Compute leading eigen pair of modularity Matrix B_hat[g].
 *  Returns the leading eigen value.
 * @param Bg - a modularity matrix.
 * @param leadEigenVec - a preallocated vector pointer, at which the leading eigenvector (corresponding to the leading eigenvalue returned) will be stored.
 */
scalar leading_eigenpair(modMat *Bg, vector *leadEigenVec);


#endif /* LEADING_EIGENPAIR_H_ */