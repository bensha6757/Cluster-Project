#ifndef LEADING_EIGENPAIR_H_
#define LEADING_EIGENPAIR_H_


#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Types.h"

#define IS_POSITIVE(X) ((X) > 0.00001)


/**
 * Compute dot product of two vectors of length d.
 */
double dot_prod(vector v, vector u, num d);

/** Compute leading eigen pair of modularity Matrix B_hat[g], and store results in referenced pointers.
 * @param Bg - a modularity matrix.
 * @param leadEigenVec - the address in memory to store the leading eigenvector (vector) in.
 * @param leadEigenVal - the address in memory to store the leading eigenvalue (scalar) in.
 */
void leading_eigenpair(modMat *Bg, vector *leadEigenVec, scalar *leadEigenVal);


#endif /* LEADING_EIGENPAIR_H_ */