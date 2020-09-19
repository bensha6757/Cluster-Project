#ifndef Leading_eigenpair_H_
#define Leading_eigenpair_H_


#include "IO_Mem_Errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Types.h"

/**********************************************************************************
 * A module in charge of supplying the leading eigenpair for Algorithm 2 purposes.*
 * ********************************************************************************/

 /*There are about O(N^2) Power iterations for a matrix of size N */
#define ITER_LIMIT(N) (0.5*(N)*(N)+2000*(N)+50000)


/**
 * Returns a dot product (scalar) of two vectors of dimension d.
 */
double dot_prod(vector v, vector u, num d);

/** Compute leading eigen pair of modularity Matrix B_hat[g].
 *  Returns the leading eigen value.
 * @param Bg - a modularity matrix.
 * @param leadEigenVec - a preallocated vector pointer, at which the leading eigenvector (corresponding to the leading eigenvalue returned) will be stored.
 */
scalar Leading_eigenpair(modMat *Bg, vector *leadEigenVec);


#endif /* Leading_eigenpair_H_ */