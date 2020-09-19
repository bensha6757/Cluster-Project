#ifndef DIVIDE_INTO_TWO_H_
#define DIVIDE_INTO_TWO_H_

#include "IO_Mem_Errors.h"
#include "modMat.h"
#include "Leading_eigenpair.h"
#include "Types.h"
#include "Flag_Set.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/**********************************************************************************
 * A module in charge of dividing a network into two sub-networks. 				  *
 **********************************************************************************/

typedef enum div_res_t {
	GROUP_INDIVISIBLE,
	GROUP_DIVIDED
} DIV_RESULT;

/** Divides a sub-network of vertices indexed in Subgroup g into a 2-partition.
 * 	Applies several computational methods based on modularity matrix of the parent network.
 *	Stores the 2-partition found in given arguments g1, g2, sizeG1, sizeG2, s.t:
 * 	union(g1,g2).isEqual(g) && sizeG1 + sizeG2 == sizeG
 * 	
 *	Either sizeG1 or sizeG2 might be zero at the end of the function. resulting effctively in:
 * 	g1.isEqual(g) ^ g2.isEqual(g)
 * 
 * 	Implements Algorithm 2 of the project.
 * 
 * 	The function also returns an enum indicating whether a spectral division is feasable.
 * 
 * 	@param B - Modularity matrix of the parent network.
 * 	@param g - An ordered subset of 0,1,...,n where n==|V| of network.
 * 	@param sizeG - size of Subgroup g (number of indices).
 * 	@param g1 - An address for storing output Subgroup g1.
 * 	@param sizeG1 - An address for storing the size of output Subgroup g1.
 *  @param g2 - An address for storing output Subgroup g2.
 * 	@param sizeG2 - An address for storing the size of output Subgroup g2.
 */
DIV_RESULT divide_into_two(modMat *B, Subgroup g, num sizeG, Subgroup *g1, Subgroup *g2, num *sizeG1, num *sizeG2);

#endif /* DIV_INTO_TWO_H_ */