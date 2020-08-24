/*
 * leading_eigenpair.h
 *
 *  Created on: 24 ����� 2020
 *      Author: ��
 */

#ifndef LEADING_EIGENPAIR_H_
#define LEADING_EIGENPAIR_H_


#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>

#define IS_POSITIVE(X) ((X) > 0.00001)

typedef enum message_t {
	GROUP_INDIVISIBLE,
	GROUP_DIVIDED
} MESSAGE;

typedef double* vector;
typedef unsigned int boolean;


/*
 * Compute leading eigen pair of modularity Matrix B_hat[g]
 */
void leading_eigenpair(modMat *B, modMat *Bg, vector leadEigenVec, double* leadEigenVal);


#endif /* LEADING_EIGENPAIR_H_ */
