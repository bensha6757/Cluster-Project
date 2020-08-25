/*
 * divide_into_two.h
 *
 *  Created on: 24 ����� 2020
 *      Author: ��
 */

#ifndef DIVIDE_INTO_TWO_H_
#define DIVIDE_INTO_TWO_H_

#include "io_mem_errors.h"
#include "modMat.h"
#include <stdio.h>
#include <stdlib.h>

#define INITIAL_Q -1

typedef enum divres_t {
	GROUP_INDIVISIBLE,
	GROUP_DIVIDED
} DIV_RESULT;

DIV_RESULT div_into_two(modMat *B,Subgroup g, size_t sizeG, Subgroup **g1, Subgroup **g2, size_t *sizeG1, size_t *sizeG2);

#endif /* DIV_INTO_TWO_H_ */
